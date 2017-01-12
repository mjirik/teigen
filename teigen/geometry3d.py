import numpy as np
from scipy.spatial.distance import cdist


def translate(point, vector, length=None):
    vector = np.asarray(vector)
    if length is not None:
        vector = length * vector / np.linalg.norm(vector)
    return (np.asarray(point) + vector).tolist()



def closest_node_2d(node, nodes, return_more=False):
    """

    :param node:
    :param nodes:
    :param return_more: return closest_node, id, dist
    :return:
    """
    dst = cdist([node], nodes)
    id = dst.argmin()
    if return_more:
        return nodes[id], id, dst.min()
    return nodes[id]

def closest_node(*args, **kwargs):
    dist_2 = node_to_spheres_square_dist(*args, **kwargs)
    return np.argmin(dist_2)

def closest_node_square_dist(*args, **kwargs):
    dist_2 = node_to_spheres_square_dist(*args, **kwargs)
    return np.min(dist_2)

def node_to_spheres_square_dist(node, nodes, nodes_radius=None):
    """
    return point distance to spheres surface

    :param node: one point
    :param nodes: center
    :param nodes_radius:
    :return:
    """
    nodes = np.asarray(nodes)
    dist_2 = np.sum((nodes - node)**2, axis=1)
    if nodes_radius is not None:
        dist_2 = dist_2 - np.asarray(nodes_radius)**2
    return dist_2

def get_points_in_line_segment(nodeA, nodeB, step): #, radius, cylinder_id):
    nodeA = np.asarray(nodeA)
    nodeB = np.asarray(nodeB)
    nodes = []
    nodes.append(nodeA)
    vector = (nodeA - nodeB).tolist()
    dist = np.linalg.norm(vector)
    while dist > step:
        nodeA = translate(nodeA, vector, -step)
        nodes.append(nodeA)
        vector = (nodeA - nodeB).tolist()
        dist = np.linalg.norm(vector)
    nodes.append(nodeB)
    return nodes


def circle(center, perp_vect, radius, element_number=10):
    """
    Function computed the circle points. No drawing.
    perp_vect is vector perpendicular to plane of circle
    """
    # tl = [0, 0.2, 0.4, 0.6, 0.8]
    tl = np.linspace(0, 1, element_number)

    # vector form center to edge of circle
    # u is a unit vector from the centre of the circle to any point on the
    # circumference

    # normalized perpendicular vector
    n = perp_vect / np.linalg.norm(perp_vect)

    # normalized vector from the centre to point on the circumference
    u = perpendicular_vector(n)
    u = u / np.linalg.norm(u)

    pts = []

    for t in tl:
        # u = np.array([0, 1, 0])
        # n = np.array([1, 0, 0])
        pt = radius * np.cos(t * 2 * np.pi) * u +\
            radius * np.sin(t * 2 * np.pi) * np.cross(u, n) +\
            center

        pt = pt.tolist()
        pts.append(pt)

    return pts


def perpendicular_vector(v):
    r""" Finds an arbitrary perpendicular vector to *v*."""
    if v[1] == 0 and v[2] == 0:
        if v[0] == 0:
            raise ValueError('zero vector')
        else:
            return np.cross(v, [0, 1, 0])
    return np.cross(v, [1, 0, 0])


def cylinder_circles(nodeA, nodeB, radius, element_number=10):
    """
    Return list of two circles with defined parameters.
    """

    vector = (np.array(nodeA) - np.array(nodeB)).tolist()
    ptsA = circle(nodeA, vector, radius, element_number)
    ptsB = circle(nodeB, vector, radius, element_number)

    return ptsA, ptsB

def plane_fit(points):
    """
    p, n = plane_fit(points)

    Given an array, points, of shape (d,...)
    representing points in d-dimensional space,
    fit an d-dimensional plane to the points.
    Return a point, p, on the plane (the point-cloud centroid),
    and the normal, n.
    """
    import numpy as np
    from numpy.linalg import svd
    points = np.reshape(points, (np.shape(points)[0], -1)) # Collapse trialing dimensions
    assert points.shape[0] <= points.shape[1], "There are only {} points in {} dimensions.".format(points.shape[1], points.shape[0])
    ctr = points.mean(axis=1)
    x = points - ctr[:,np.newaxis]
    M = np.dot(x, x.T) # Could also use np.cov(x) here.
    return ctr, svd(M)[0][:,-1]

def is_in_area(pt, areasize_px, radius=None):
    """
    check if point is in area with considering eventual maximum radius
    :param node:
    :param radius:
    :return:
    """
    if areasize_px is None:
        return True

    node = np.asarray(pt)
    if radius is None:
        radius = 0

    if np.all(node > (0 + radius)) and np.all(node < (areasize_px - radius)):
        return  True
    else:
        return False

def is_cylinder_in_area(pt1, pt2, radius, areasize):
    return is_in_area(pt1, areasize, radius) and is_in_area(pt2, areasize,  radius)

def cylinder_collision(
        pt1_mm,
        pt2_mm,
        radius_mm,
        other_points,
        other_points_radiuses=None,
        areasize_mm=None,
        DIST_MAX_RADIUS_MULTIPLICATOR=1.5, # higher than sqrt(2)
        OVERLAPS_ALOWED=False
):

    step = 2 * radius_mm

    if pt1_mm is not None and is_cylinder_in_area(pt1_mm, pt2_mm, radius_mm, areasize_mm):
        line_nodes = get_points_in_line_segment(pt1_mm, pt2_mm, step)
        if OVERLAPS_ALOWED:
            return False, line_nodes

        if len(other_points) == 0:
            return False, line_nodes
        else:
            safe_dist2 = (radius_mm * DIST_MAX_RADIUS_MULTIPLICATOR) ** 2
            for node in line_nodes:
                dist_closest = closest_node_square_dist(node, other_points, other_points_radiuses)
                if dist_closest < safe_dist2:
                    return True, []
            return False, line_nodes
    return True, []
