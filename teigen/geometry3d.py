#! /usr/bin/env python
# -*- coding: utf-8 -*-

import logging
logger = logging.getLogger(__name__)

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

def n_closest_nodes(*args, **kwargs):
    """

    :param args:
    :param kwargs:
    :return:
    """
    if "n" in kwargs.keys():
        n = kwargs.pop("n")
    dist_2 = node_to_spheres_dist(*args, **kwargs)
    indexes = np.argsort(dist_2)[:n]
    distances = dist_2[indexes]

    return indexes, distances

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
    pts, step = get_points_in_line_segment(pt1, pt2, step=radius*relative_step, limit_step_number=1000)
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


def closest_distance_between_lines(a0,a1,b0,b1,clampAll=False,clampA0=False,clampA1=False,clampB0=False,clampB1=False):

    ''' Given two lines defined by numpy.array pairs (a0,a1,b0,b1)
        Return the closest points on each segment and their distance,

        Based on Fnord implementation. See:
        http://stackoverflow.com/questions/2824478/shortest-distance-between-two-line-segments

    '''

    # If clampAll=True, set all clamps to True
    if clampAll:
        clampA0=True
        clampA1=True
        clampB0=True
        clampB1=True

    a0 = np.asarray(a0)
    a1 = np.asarray(a1)
    b0 = np.asarray(b0)
    b1 = np.asarray(b1)

    # Calculate denomitator
    A = a1 - a0
    B = b1 - b0
    magA = np.linalg.norm(A)
    magB = np.linalg.norm(B)

    _A = A / magA
    _B = B / magB

    cross = np.cross(_A, _B);
    denom = np.linalg.norm(cross)**2


    # If lines are parallel (denom=0) test if lines overlap.
    # If they don't overlap then there is a closest point solution.
    # If they do overlap, there are infinite closest positions, but there is a closest distance
    if not denom:
        d0 = np.dot(_A,(b0-a0))

        # Overlap only possible with clamping
        if clampA0 or clampA1 or clampB0 or clampB1:
            d1 = np.dot(_A,(b1-a0))

            # Is segment B before A?
            if d0 <= 0 >= d1:
                if clampA0 and clampB1:
                    if np.absolute(d0) < np.absolute(d1):
                        return a0,b0,np.linalg.norm(a0-b0)
                    return a0,b1,np.linalg.norm(a0-b1)


            # Is segment B after A?
            elif d0 >= magA <= d1:
                if clampA1 and clampB0:
                    if np.absolute(d0) < np.absolute(d1):
                        return a1,b0,np.linalg.norm(a1-b0)
                    return a1,b1,np.linalg.norm(a1-b1)


        # Segments overlap, return distance between parallel segments
        return None,None,np.linalg.norm(((d0*_A)+a0)-b0)



    # Lines criss-cross: Calculate the projected closest points
    t = (b0 - a0);
    detA = np.linalg.det([t, _B, cross])
    detB = np.linalg.det([t, _A, cross])

    t0 = detA/denom;
    t1 = detB/denom;

    pA = a0 + (_A * t0) # Projected closest point on segment A
    pB = b0 + (_B * t1) # Projected closest point on segment B


    # Clamp projections
    if clampA0 or clampA1 or clampB0 or clampB1:
        if clampA0 and t0 < 0:
            pA = a0
        elif clampA1 and t0 > magA:
            pA = a1

        if clampB0 and t1 < 0:
            pB = b0
        elif clampB1 and t1 > magB:
            pB = b1

        # Clamp projection A
        if (clampA0 and t0 < 0) or (clampA1 and t0 > magA):
            dot = np.dot(_B,(pA-b0))
            if clampB0 and dot < 0:
                dot = 0
            elif clampB1 and dot > magB:
                dot = magB
            pB = b0 + (_B * dot)

        # Clamp projection B
        if (clampB0 and t1 < 0) or (clampB1 and t1 > magB):
            dot = np.dot(_A,(pB-a0))
            if clampA0 and dot < 0:
                dot = 0
            elif clampA1 and dot > magA:
                dot = magA
            pA = a0 + (_A * dot)


    return pA,pB,np.linalg.norm(pA-pB)

def bbox_collision(bbox1, bbox2):
    """
    detects collision betwen two boundingboxes.

    :param bbox1: [[minX, maxX], [minY, maxY] ... [...]]
    :param bbox2: [[minX, maxX], [minY, maxY] ... [...]]
    :return:
    """

    bbox1 = np.asarray(bbox1)
    bbox2 = np.asarray(bbox2)

    max1 = np.max(bbox1, axis=1)
    min1 = np.min(bbox1, axis=1)

    max2 = np.max(bbox2, axis=1)
    min2 = np.min(bbox2, axis=1)

    out = (min1 <= max2) & (max1 >= min2)
    return np.all(out)

def get_bbox(points, margin=0):
    """
    Get bounding box based on points. Margin can be added
    :param points:
    :param margin:
    :return:
    """
    points = np.asarray(points)
    ptsplus = points + margin
    ptsminus = points - margin

    points = np.concatenate((ptsplus, ptsminus))
    bbox = np.zeros([points.shape[1], 2])
    bbox[:, 0] = np.min(points, axis=0)
    bbox[:, 1] = np.max(points, axis=0)

    return bbox

def get_bbox_corners(bbox):
    import itertools

    pts = []
    for prod in itertools.product([0, 1], repeat=len(bbox)):
        pt = np.zeros([len(prod)])
        for axi in range(len(prod)):
            pt[axi] = bbox[axi][prod[axi]]
        pts.append(pt)

    return pts



def cylinder_collision_detection(pointA1, pointA2, radiusA, pointB1, pointB2, radiusB, bboxA=None, bboxB=None):
    """
    Detects the possibility of cylinder intersection.
gitggggglkjlkjaasdfasdf
    :param pointA1:
    :param pointA2:
    :param radiusA:
    :param pointB1:
    :param pointB2:
    :param radiusB:
    :param bboxA:
    :param bboxB:
    :return:
    """

    if bboxA is None:
        bboxA = get_bbox([pointA1, pointA2], margin=radiusA)
    if bboxB is None:
            bboxB = get_bbox([pointB1, pointB2], margin=radiusB)

def point_and_plane(plane_point, plane_orientation, points=None, xyz=None):
    """
    return 0 if point is in plane
    :param plane_point:
    :param plane_orientation:
    :param points: [[x1, y1, z1], [x2, y2, z2], ... ]
    :param xyz: [x,y,z] or more points [[x1, x2, ...], [y1, y2, ...], [z1, z2, ...]]
    :return:
    """
    vector = plane_orientation
    vector = vector / np.linalg.norm(vector)
    a = vector[0]
    b = vector[1]
    c = vector[2]

    d = -a * plane_point[0] - b * plane_point[1] - c * plane_point[2]

    if xyz is not None:
        xyz = np.asarray(xyz)
    elif points is not None:
        points = np.asarray(points)
        xyz = points.T
    else:
        logger.error("points or xyz must be declared")


    x, y, z = xyz
    z_out = (a * x + b * y + c*z + d) / (a**2 + b**2 +c**2)**0.5

    return z_out


class GeometricObject():
    def __init__(self, bbox=None):
        self.bbox = bbox

    def bbox_collision(self, bbox):
        return bbox_collision(self.bbox, bbox)

class CylinderObject(GeometricObject):

    def __init__(self, point1, point2, radius):

        bbox = get_bbox([point1, point2], margin=radius)
        GeometricObject.__init__(self, bbox=bbox)
        self.point1 = point1
        self.point2 = point2
        self.radius = radius

    def _separable_by_bases(self, obj):
        sep1 = self._separable_by_one_bbox_and_base(obj.bbox, self.point1, self.point2)
        sep2 = self._separable_by_one_bbox_and_base(obj.bbox, self.point2, self.point1)
        sep3 = self._separable_by_one_bbox_and_base(self.bbox, obj.point2, obj.point1)
        sep4 = self._separable_by_one_bbox_and_base(self.bbox, obj.point2, obj.point1)
        return sep1 | sep2 | sep3 | sep4

    def _separable_by_one_bbox_and_base(self, bbox, base_point, other_point):
        vector = np.asarray(other_point) - np.asarray(base_point)
        points = get_bbox_corners(bbox)
        position = point_and_plane(base_point, vector, points)
        return np.all(position < 0)

    def _separable_by_dist(self, obj):

        dist = closest_distance_between_lines(obj.point1, obj.point2, self.point1, self.point2, clampAll=True)
        return dist > (obj.radius + self.radius)

    def _separable_by_bbox(self, obj):
        return not self.bbox_collision(obj.bbox)

    def collision(self, obj):
        if self._separable_by_bbox(obj):
            # are separable by bbox
            return False
        else:
            # TODO Implement type check
            # if type(obj) == CylinderObject:
            if True:
                if self._separable_by_dist(obj):
                    return False
                elif self._separable_by_bases(obj):
                    return False
        return True


class CollisionModel():
    pass

class CollisionSpheresModel():

    def __init__(self, areasize=None):
        self.collision_alowed = False
        self._cylinder_nodes = []
        self._cylinder_nodes_radiuses = []
        self._cylinder_end_nodes = []
        self._cylinder_end_nodes_radiuses = []
        if areasize is not None:
            areasize = np.asarray(areasize)
        self.areasize = areasize
        self.object_number = 0


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
            self._cylinder_end_nodes.append(pt1)
            self._cylinder_end_nodes.append(pt2)
            self._cylinder_end_nodes_radiuses.append(radius)
            self._cylinder_end_nodes_radiuses.append(radius)
            self.object_number += 1

        return collision

    def get_random_point(self, radius=None):
        if radius is not None:
            pt1 = (np.random.random([3]) * (self.areasize - (2 * radius))) + radius
        else:
            pt1 = np.random.random([3]) * self.areasize

        return pt1

    def get_node_number(self):
        return len(self._cylinder_nodes)

    def n_closest_points(self, node, n):
        indexes, distances = n_closest_nodes(node=node, n=n, nodes=self._cylinder_nodes, nodes_radius=self._cylinder_nodes_radiuses)
        nodes = np.asarray(self._cylinder_nodes)[indexes]
        return nodes, indexes, distances

    def n_closest_end_points(self, node, n):
        indexes, distances = n_closest_nodes(node=node, n=n, nodes=self._cylinder_end_nodes, nodes_radius=self._cylinder_end_nodes_radiuses)
        nodes = np.asarray(self._cylinder_end_nodes)[indexes]
        return nodes, indexes, distances
