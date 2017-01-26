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
    dist_2 = node_to_spheres_dist(*args, **kwargs)
    return np.argmin(dist_2)

def closest_node_dist(*args, **kwargs):
    dist_2 = node_to_spheres_dist(*args, **kwargs)
    min_dst_2 = np.min(dist_2)
    return min_dst_2

def node_to_spheres_dist(node, nodes, nodes_radius=None, return_square=False):
    """
    return point distance to spheres surface

    :param node: one point
    :param nodes: center
    :param nodes_radius:
    :param return_square: faster but no sqrt is performed
    :return:
    """
    vectors = np.asarray(nodes) - np.asarray(node)
    # dist = np.sum(vectors**2, axis=1)**0.5
    # dist = np.sum((np.asarray(nodes) - node)**2, axis=1)**0.5
    dist = np.linalg.norm(vectors, axis=1)

    nodes = np.asarray(nodes)
    if nodes_radius is not None:
        dist = dist - np.asarray(nodes_radius)
    return dist

def get_spheres_bounding_cylinder(pt1, pt2, radius): #, relative_step=0.5):
    # step to raidus
    relative_step = 0.5
    safety = 1.00001
    # sphere_radius_ratio = ((1 + relative_step**2)**0.5) * safety
    sphere_radius_ratio = 1.118034

    # relative_step = 1.0
    # constant higher than sqrt(2)
    # sphere_radius_ratio = 1.414214
    pts, step = get_points_in_line_segment(pt1, pt2, step=radius*relative_step, limit_step_number=100)
    if radius == step:
        radiuses = [radius * sphere_radius_ratio] * len(pts)
    else:
        radiuses = [(step**2 + radius**2)**0.5 * safety] * len(pts)
    return pts, radiuses

def get_points_closer(nodeA, nodeB, delta=None, relative_length=None): #, radius, cylinder_id):
    vector = (np.asarray(nodeA) - np.asarray(nodeB)).tolist()
    length = np.linalg.norm(vector)

    if relative_length is not None:
        delta = 0.5 * (length - (length * relative_length))
    else:
        delta = 0.5 * delta


    if length < 2*delta:
        return None, None

    # mov circles to center of cylinder by size of radius because of joint
    nodeA = translate(nodeA, vector,
                         -delta) # * self.endDistMultiplicator)
    nodeB = translate(nodeB, vector,
                         delta) #  * self.endDistMultiplicator)
    nodeA = np.asarray(nodeA)
    nodeB = np.asarray(nodeB)
    return nodeA, nodeB

def get_points_in_line_segment(nodeA, nodeB, step, limit_step_number=None): #, radius, cylinder_id):
    nodeA = np.asarray(nodeA)
    nodeB = np.asarray(nodeB)
    nodes = []
    nodes.append(nodeA)
    vector = (nodeA - nodeB).tolist()
    dist = np.linalg.norm(vector)
    if limit_step_number is not None:
        if dist/step > limit_step_number:
            step = dist * 1. / limit_step_number
    while dist > step:
        nodeA = translate(nodeA, vector, -step)
        nodes.append(nodeA)
        vector = (nodeA - nodeB).tolist()
        dist = np.linalg.norm(vector)
    nodes.append(nodeB)

    if limit_step_number is not None:
        return nodes, step
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

def is_in_area(pt, areasize, radius=None):
    """
    check if point is in area with considering eventual maximum radius
    :param node:
    :param radius:
    :return:
    """
    if areasize is None:
        return True

    node = np.asarray(pt)
    if radius is None:
        radius = 0

    if np.all(node > (0 + radius)) and np.all(node < (areasize - radius)):
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
        # DIST_MAX_RADIUS_MULTIPLICATOR=1.414214, # higher than sqrt(2)
        COLLISION_ALOWED=False
):

    if pt1_mm is not None and is_cylinder_in_area(pt1_mm, pt2_mm, radius_mm, areasize_mm):
        # line_nodes = get_points_in_line_segment(pt1_mm, pt2_mm, step)
        line_nodes, nodes_radiuses = get_spheres_bounding_cylinder(pt1_mm, pt2_mm, radius=radius_mm)
        if COLLISION_ALOWED:
            return False, line_nodes, nodes_radiuses

        if len(other_points) == 0:
            return False, line_nodes, nodes_radiuses
        else:
            # safe_dist2 = radius_mm * DIST_MAX_RADIUS_MULTIPLICATOR
            for node, safe_dist in zip(line_nodes, nodes_radiuses):
                dist_closest = closest_node_dist(node, other_points, other_points_radiuses)
                if dist_closest < safe_dist:
                    return True, [], []
            return False, line_nodes, nodes_radiuses
    return True, [], []


class CollisionBoundaryModel():

    def __init__(self, areasize=None):
        self.collision_alowed = False
        self._cylinder_nodes = []
        self._cylinder_nodes_radiuses = []
        if areasize is not None:
            areasize = np.asarray(areasize)
        self.areasize = areasize


    def is_in_area(self, node, radius=None):
        """
        check if point is in area with considering eventual maximum radius
        :param node:
        :param radius:
        :return:
        """
        node = np.asarray(node)
        if radius is None:
            radius = self.radius_maximum
        return is_in_area(node, self.areasize, radius=radius)

    def add_cylinder_if_no_collision(self, pt1, pt2, radius,
                            # COLLISION_RADIUS=1.5 # higher then sqrt(2)
                            ):
        # TODO use geometry3.check_collision_along_line
        collision, new_nodes, nodes_radiuses = cylinder_collision(
            pt1,
            pt2,
            radius,
            other_points=self._cylinder_nodes,
            other_points_radiuses=self._cylinder_nodes_radiuses,
            areasize_mm=self.areasize,
            # DIST_MAX_RADIUS_MULTIPLICATOR=self.DIST_MAX_RADIUS_MULTIPLICATOR,
            COLLISION_ALOWED=self.collision_alowed,
        )

        if not collision:
            self._cylinder_nodes.extend(new_nodes)
            self._cylinder_nodes_radiuses.extend(nodes_radiuses)

        return collision

    def get_random_point(self, radius=None):
        if radius is not None:
            pt1 = (np.random.random([3]) * (self.areasize - (2 * radius))) + radius
        else:
            pt1 = np.random.random([3]) * self.areasize

        return pt1
