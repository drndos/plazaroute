import itertools
import logging
from typing import List

import rtree
from shapely.geometry import CAP_STYLE
from shapely.geometry import Point, MultiPolygon, Polygon, LineString, box
from shapely.ops import unary_union

from plaza_preprocessing import configuration
from plaza_preprocessing.importer.osmholder import OSMHolder
from plaza_preprocessing.optimizer import utils
from plaza_preprocessing.optimizer.graphprocessor.graphprocessor import \
    GraphProcessor

logger = logging.getLogger('plaza_preprocessing.optimizer')


def preprocess_plazas(osm_holder: OSMHolder, process_strategy: GraphProcessor, shortest_path_strategy, config: dict):
    """ preprocess all plazas from osm_importer """
    logger.info(f"Start processing {len(osm_holder.plazas)} plazas")
    plaza_processor = PlazaPreprocessor(
        osm_holder, process_strategy, shortest_path_strategy, config)
    processed_plazas = plaza_processor.process_plazas()

    logger.info(f"Finished processing {len(processed_plazas)} plazas (rest were discarded)")
    return processed_plazas


class PlazaPreprocessor:

    def __init__(self, osm_holder: OSMHolder, graph_processor: GraphProcessor,
                 shortest_path_strategy, config):
        self.plazas = osm_holder.plazas
        self.lines = osm_holder.lines
        self.buildings = osm_holder.buildings
        self.points = osm_holder.points
        self.graph_processor = graph_processor
        self.shortest_path_strategy = shortest_path_strategy
        self.config = config
        self.no_plaza_intersections = {}

        self._create_spatial_indices()

    def process_plazas(self):
        """ process all plazas in the osm holder"""
        processed_plazas = []
        merged_plazas = self._merge_plazas(self.plazas)
        for key, plaza in merged_plazas.items():
            logger.info(f"Processing plaza {plaza['osm_id']}")
            processed_plaza = self._process_plaza(plaza)
            if processed_plaza is not None:
                processed_plazas.append(processed_plaza)

        return processed_plazas

    def _merge_plazas(self, plazas):
        merged_plazas = {}
        intersecting_plazas = {}
        osm_id = 0

        for plaza in plazas:
            merged_plazas[plaza['osm_id']] = plaza

        for plaza in plazas:
            intersecting_plazas[plaza['osm_id']] = plaza

        while True:
            intersections = self._find_intersections(intersecting_plazas)
            intersecting_plazas = {}
            if len(intersections) > 0:
                for intersection in intersections:
                    if intersection['left']['osm_id'] in merged_plazas or intersection['right']['osm_id'] in merged_plazas:
                        osm_id = osm_id + 1
                        logger.info(f"Merging plazas {intersection['left']['osm_id']} and {intersection['right']['osm_id']} into {osm_id}"),
                        to_merge = [intersection['left']['geometry'], intersection['right']['geometry']]
                        merged_plaza_geometry = unary_union(to_merge)
                        merged_plaza = {
                            'osm_id': osm_id,
                            'geometry': merged_plaza_geometry
                        }
                        merged_plazas.pop(intersection['left']['osm_id'], None)
                        merged_plazas.pop(intersection['right']['osm_id'], None)
                        merged_plazas[osm_id] = merged_plaza
                        intersecting_plazas[osm_id] = merged_plaza
            else:
                break

        return merged_plazas

    def _find_intersections(self, plazas: dict):
        intersections = []
        combinations = list(itertools.combinations(list(plazas.values()), 2))
        for comb in combinations:
            logger.debug(f"Looking for intersection {comb[0]} and {comb[1]}")
            if comb[0]['geometry'].intersects(comb[1]['geometry']):
                intersections.append({'left': comb[0], 'right': comb[1]})

        return intersections

    def _create_spatial_indices(self):
        """ create spatial indices for lines, buildings and points"""
        logger.info("Creating spatial index for geometries")
        line_geometries = [line['geometry'] for line in self.lines]
        self.line_index = self._create_spatial_index(line_geometries)
        self.building_index = self._create_spatial_index(self.buildings)
        self.point_index = self._create_spatial_index(self.points)

    def _process_plaza(self, plaza):
        """ process a single plaza """

        intersecting_lines = self._find_intersecting_lines(plaza['geometry'])

        plaza_geom_without_obstacles = self._calc_obstacle_geometry(
            plaza, intersecting_lines, buffer_m=self.config['obstacle-buffer'])

        if not plaza_geom_without_obstacles:
            logger.debug(f"Discarding Plaza {plaza['osm_id']}: completely obstructed by obstacles")
            return None

        entry_points = self._calc_entry_points(plaza_geom_without_obstacles, intersecting_lines, lookup_buffer_m=self.config['entry-point-lookup-buffer'])

        if len(entry_points) < 2:
            logger.debug(f"Discarding Plaza {plaza['osm_id']} - it has fewer than 2 entry points")
            return None

        entry_lines = self._map_entry_lines(intersecting_lines, entry_points)

        graph_edges = self._get_graph_edges(entry_points, plaza['geometry'], plaza_geom_without_obstacles)

        if not graph_edges:
            logger.debug(f"Discarding Plaza {plaza['osm_id']}: no graph could be constructed")
            return None

        plaza['geometry'] = plaza_geom_without_obstacles
        plaza['entry_points'] = entry_points
        plaza['entry_lines'] = entry_lines
        plaza['graph_edges'] = graph_edges

        return plaza

    def _get_graph_edges(self, entry_points: List[Point], plaza_geom: Polygon,
                         plaza_geom_without_obstacles: Polygon) -> List[LineString]:
        """ create graph with shortest paths between entry points """

        """
         .*\(([^ ]+) ([^)]+)\)
         L.marker([$2, $1]).addTo(map);
        """
        graph_edges = self.graph_processor.create_graph_edges(plaza_geom_without_obstacles, entry_points)

        return graph_edges

    def _calc_entry_points(self, plaza_geometry, intersecting_lines, lookup_buffer_m):
        """
        calculate points where lines intersect with the outer ring of the plaza
        """
        intersection_coords = set()
        for line in intersecting_lines:
            line_geom = line['geometry']
            intersection = line_geom.intersection(plaza_geometry)
            line_intersection_coords = utils.unpack_geometry_coordinates(intersection)
            intersection_coords = intersection_coords.union(line_intersection_coords)

        intersection_points = list(map(Point, intersection_coords))

        return intersection_points

    def _map_entry_lines(self, intersecting_lines, entry_points):
        """ map entry lines to entry points """
        entry_lines = []
        for line in intersecting_lines:
            matching_entry_points = list(filter(
                lambda p: (p.x, p.y) in line['geometry'].coords, entry_points))
            if matching_entry_points:
                entry_lines.append({
                    'way_id': line['id'],
                    'entry_points': matching_entry_points
                })
        return entry_lines

    def _find_intersecting_lines(self, plaza_geometry):
        """ return every line that intersects with the plaza """
        intersecting_lines = []
        potential_matches = self._search_index(
            self.line_index, plaza_geometry.bounds, self.lines)
        for line in potential_matches:
            if plaza_geometry.intersects(line['geometry']):
                intersecting_lines.append(line)

        return intersecting_lines

    def _calc_obstacle_geometry(self, plaza, intersecting_lines, buffer_m):
        """ cuts out holes for obstacles on the plaza geometry """
        intersecting_buildings = self._find_intersecting_buildings(plaza['geometry'])

        geometry_without_buildings = plaza['geometry']
        for building in intersecting_buildings:
            geometry_without_buildings = geometry_without_buildings.difference(building)

        points_on_plaza = self._get_points_inside_plaza(plaza['geometry'])
        point_obstacles = list(
            map(lambda p: self._create_point_obstacle(p, buffer_m), points_on_plaza))

        barrier_obstacles = self._create_barrier_obstacles(intersecting_lines, self.config['obstacle-buffer'] / 2)

        geometry_without_obstacles = geometry_without_buildings
        for point_obstacle in point_obstacles:
            geometry_without_obstacles = geometry_without_obstacles.difference(point_obstacle)
        for barrier_obstacle in barrier_obstacles:
            geometry_without_obstacles = geometry_without_obstacles.difference(barrier_obstacle)

        if isinstance(geometry_without_obstacles, MultiPolygon):
            logger.debug(
                f"Plaza {plaza['osm_id']}: Multipolygon after cut out, discarding smaller polygon")
            # take the largest of the polygons
            largest_geometry_without_obstacles = max(
                geometry_without_obstacles, key=lambda p: p.area)
            # if cut out is less than 5% the area of the original, it's discarded
            if largest_geometry_without_obstacles.area < plaza['geometry'].area * 0.05:
                return None
            return largest_geometry_without_obstacles

        return geometry_without_obstacles

    def _find_intersecting_buildings(self, plaza_geometry):
        """ finds all buildings on the plaza that have not been cut out"""
        potential_matches = self._search_index(
            self.building_index, plaza_geometry.bounds, self.buildings)
        return list(filter(plaza_geometry.intersects, potential_matches))

    def _get_points_inside_plaza(self, plaza_geometry):
        """ finds all points that are on the plaza geometry """
        potential_matches = self._search_index(
            self.point_index, plaza_geometry.bounds, self.points)
        return list(filter(plaza_geometry.intersects, potential_matches))

    def _create_spatial_index(self, geometries):
        """ create rtree index for fast intersection checking """
        logger.debug(f"creating spatial index for {len(geometries)} geometries")
        idx = rtree.index.Index()
        for i, geometry in enumerate(geometries):
            idx.insert(i, geometry.bounds)
        return idx

    def _search_index(self, index, bounds, geometries):
        """
        search rtree index and return geometries that potentially
        intersect with the bounds
        """
        potential_matches_indices = index.intersection(bounds)
        return map(lambda i: geometries[i], potential_matches_indices)

    def _create_point_obstacle(self, point, buffer_m):
        """ create a polygon around a point with a buffer in meters """
        buffer_deg = utils.meters_to_degrees(buffer_m)
        min_x = point.x - buffer_deg
        min_y = point.y - buffer_deg
        max_x = point.x + buffer_deg
        max_y = point.y + buffer_deg
        return box(min_x, min_y, max_x, max_y)

    def _create_barrier_obstacles(self, intersecting_lines, buffer_m):
        """ returns geometries for line obstacles, e.g. barriers"""
        tag_filter = self.config['tag-filter']['barrier']
        buffer_distance = utils.meters_to_degrees(buffer_m)
        barrier_obstacles = filter(lambda line: configuration.filter_tags(line['tags'], tag_filter), intersecting_lines)
        buffered_obstacles = map(
            lambda l: l['geometry'].buffer(buffer_distance, cap_style=CAP_STYLE.flat), barrier_obstacles)
        return buffered_obstacles
