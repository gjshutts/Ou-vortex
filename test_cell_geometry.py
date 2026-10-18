"""Known geometry, conservation of coordinate area, and envelope checks."""

import math
import unittest
from cell_geometry import (rectangle, clip_half_plane, signed_area,
                           cell_polygons, plane_value)


class CellGeometryTests(unittest.TestCase):
    def setUp(self):
        self.square = rectangle(0, 1)

    def test_diagonal_cut_has_known_vertices_and_area(self):
        triangle = clip_half_plane(self.square, -1, -1, 1)
        self.assertEqual(set(triangle), {(0, 0), (1, 0), (0, 1)})
        self.assertEqual(signed_area(triangle), 0.5)

    def test_constant_and_excluding_inequalities(self):
        for inequality in [(0, 0, -1), (1, 0, -2)]:
            self.assertEqual(clip_half_plane(self.square, *inequality), [])
        self.assertEqual(signed_area(clip_half_plane(self.square, 0, 0, 1)), 1)

    def test_boundary_contact_is_not_a_positive_area_cell(self):
        self.assertEqual(clip_half_plane(self.square, 1, 0, -1), [])
        self.assertEqual(clip_half_plane(self.square, 1, 1, -2), [])
        kept = clip_half_plane(self.square, 1, 1, 0)
        self.assertEqual(len(kept), 4)
        self.assertEqual(signed_area(kept), 1)

    def test_orientation_and_rescaling(self):
        expected = clip_half_plane(self.square, 1, -2, 0.25)
        for factor in [1e-15, 1, 1e15]:
            actual = clip_half_plane(list(reversed(self.square)),
                                     factor, -2*factor, 0.25*factor)
            self.assertGreater(signed_area(actual), 0)
            self.assertEqual(set(actual), set(expected))

    def test_partition_and_plane_dominance(self):
        planes = [(0, 0, 0), (1, 0, -0.5), (0, 1, -0.5)]
        cells = cell_polygons(planes, self.square)
        self.assertEqual([signed_area(p) for p in cells], [0.25, 0.375, 0.375])
        for i, poly in enumerate(cells):
            for point in poly:
                for plane in planes:
                    self.assertGreaterEqual(plane_value(planes[i], point) + 1e-12,
                                            plane_value(plane, point))

    def test_displaced_offsets_and_empty_cells(self):
        # P0=0 dominates P1=R-2 over the entire square.
        cells = cell_polygons([(0, 0, 0), (1, 0, -2)], self.square)
        self.assertEqual([signed_area(p) for p in cells], [1, 0])
        # Equal slopes with different offsets: only the upper plane survives.
        cells = cell_polygons([(1, 2, 0), (1, 2, 1)], self.square)
        self.assertEqual([signed_area(p) for p in cells], [0, 1])

    def test_identical_planes_do_not_duplicate_area(self):
        cells = cell_polygons([(1, 2, 0), (1, 2, 0)], self.square)
        self.assertEqual([signed_area(p) for p in cells], [1, 0])

    def test_invalid_inputs(self):
        for call in [lambda: rectangle(1, 0),
                     lambda: clip_half_plane(self.square, math.nan, 0, 0),
                     lambda: clip_half_plane(self.square, 1, 0, 0, tolerance=-1),
                     lambda: cell_polygons([(1, 2)], self.square)]:
            with self.assertRaises(ValueError):
                call()


if __name__ == "__main__":
    unittest.main(verbosity=2)
