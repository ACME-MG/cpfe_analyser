import numpy as np


def ratio(t1: np.ndarray, t2: np.ndarray) -> tuple:#[float, float]:
    """
    Args:
      t1: An array of coordinates of shape (N, 3, 2), where N is the
          number of triangles, 3 for the vertices of each triangle, and
          2 for the x and y coordinates, respectively.
      t2: Similar to t1, with shape (M, 3, 2)

    Returns a tuple of floats representing the ratio between triangle areas.
    The ratios are the minimum (rmin) and the maximum (rmax) area ratio between
    triangles in t1 and triangles in t2.
    """
    get_area = lambda x1, x2, x3, y1, y2, y3 : 0.5*abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    area1 = [get_area(*t[:,0], *t[:,1]) for t in t1]
    area2 = [get_area(*t[:,0], *t[:,1]) for t in t2]
    ratios = [a1/a2 for a1 in area1 for a2 in area2]
    return min(ratios), max(ratios)

if __name__ == "__main__":
    t1 = np.array(
        [
            [[1, 2], [2, 5], [-2, 3]],
            [[3, 2], [5, 1], [9, 2]],
            [[-3, -1], [1, 1], [7, -11]],
        ]
    )

    t2 = np.array(
        [
            [[5, 2], [7, 4], [5, -3]],
            [[-3, 2], [4, 4], [2, -13]],
        ]
    )

    rmin, rmax = ratio(t1, t2)
    assert np.isclose(rmin, 0.05217391304347826)
    assert np.isclose(rmax, 6)
