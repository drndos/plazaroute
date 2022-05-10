import logging
import osmium
from osmium._osmium import InvalidLocationError
import shapely.wkb as wkblib
from plaza_preprocessing.importer import osmholder
from plaza_preprocessing import configuration

logger = logging.getLogger('plaza_preprocessing.importer')
WKBFAB = osmium.geom.WKBFactory()

OSM_MAX_ID = 10**10


def import_osm(filename, tag_filters):
    """ imports a OSM / PBF file and returns a holder with all plazas, buildings,
    lines and points with shapely geometries """
    logger.info(f'importing {filename}')
    handler = _PlazaHandler(tag_filters)

    index_type = 'sparse_mem_array'
    handler.apply_file(filename, locations=True, idx=index_type)

    logger.debug(f'found {len(handler.plazas)} plazas')
    logger.debug(f'found {len(handler.buildings)} buildings')
    logger.debug(f'found {len(handler.lines)} lines')
    logger.debug(f'found {len(handler.points)} points')

    if handler.invalid_count > 0:
        logger.warning(f'encountered {handler.invalid_count} invalid objects (may be because of boundaries)')
    return osmholder.OSMHolder(handler.plazas, handler.buildings, handler.lines, handler.points)


class _PlazaHandler(osmium.SimpleHandler):
    def __init__(self, tag_filters):
        super().__init__()
        self.tag_filters = tag_filters
        self.plazas = []
        self.buildings = []
        self.points = []
        self.lines = []
        self.invalid_count = 0

    def node(self, node):
        if self._is_relevant_node(node):
            self._check_max_id(node.id)
            point_wkb = WKBFAB.create_point(node)
            point_geometry = wkblib.loads(point_wkb, hex=True)
            self.points.append(point_geometry)

    def way(self, way):
        if self._is_relevant_way(way):
            self._check_max_id(way.id)
            try:
                line_wkb = WKBFAB.create_linestring(way)
                line_geometry = wkblib.loads(line_wkb, hex=True)
                self.lines.append({
                    'id': way.id,
                    'geometry': line_geometry,
                    'tags': {t.k: t.v for t in way.tags}
                })
            except InvalidLocationError:
                logger.debug(f'Encountered invalid location in way {way.id}')
                self.invalid_count += 1
            except RuntimeError as ex:
                logger.debug(f'Error importing way {way.id}: {ex}')
                self.invalid_count += 1

    def area(self, area):
        if self._is_plaza(area):
            self._check_max_id(area.id)
            multipolygon_geom = self._create_multipolygon(area)
            if multipolygon_geom:
                for polygon in multipolygon_geom.geoms:
                    plaza = {
                        'osm_id': area.orig_id(),
                        'geometry': polygon
                    }
                    self.plazas.append(plaza)

        elif self._is_relevant_building(area):
            geometry = self._create_multipolygon(area)
            if geometry:
                self.buildings.append(geometry)

    def _create_multipolygon(self, area):
        try:
            building_wkb = WKBFAB.create_multipolygon(area)
            return wkblib.loads(building_wkb, hex=True)

        except InvalidLocationError:
            logger.debug(f'Encountered invalid location in area {area.id}')
            self.invalid_count += 1
            return None
        except RuntimeError as ex:
            logger.debug(f'Error importing way {area.id}: {ex}')
            self.invalid_count += 1
            return None

    def _is_relevant_node(self, node):
        return node.tags.get("level", "0") == "0" and \
            node.tags.get("layer", "0") == "0" and \
            configuration.filter_tags(node.tags, self.tag_filters['point_obstacle'])

    def _is_relevant_way(self, way):
        return "highway" in way.tags or \
            way.tags.get("railway") == "tram" or \
            configuration.filter_tags(way.tags, self.tag_filters['barrier'])

    def _is_plaza(self, area):
        return configuration.filter_tags(area.tags, self.tag_filters['plaza'])

    def _is_relevant_building(self, area):
        return "building" in area.tags \
            and area.tags.get("layer", "0") == "0"

    def _check_max_id(self, osm_id):
        if osm_id >= OSM_MAX_ID:
            logger.error(f"OSM id {osm_id} is larger than the allowed {OSM_MAX_ID}")

