from time import sleep

class Point(object):
    """docstring for Point."""

    x = 0.
    y = 0.
    z = None

    def __init__(self, x, y, z=None):
        super(Point, self).__init__()
        self.x = x
        self.y = y
        self.z = z

    def set_coords(self, x, y, z=None):
        self.x = x
        self.y = y
        self.z = z
        pass

    def get_coords(self):
        return [x, y, z]

    def get_coords2D(self):
        return [x, y]


def distance(a, b):
    dist = 0
    dist += (a.x - b.x)**2
    dist += (a.y - b.y)**2
    if a.z is not None and b.z is not None:
        dist += (a.z - b.z)**2
    return dist ** 0.5

class Region(object):
    pass


class Block(Region):
    """docstring for Block."""
    lB = Point(0., 0.)  # left Bottom point
    tR = Point(0., 0.)  # topRight point

    def __init__(self, a, b):
        super(Region, self).__init__()
        self.lB = a
        self.tR = b

    def set_Geometry(self, a, b):
        self.lB = a
        self.tR = b
        pass

    def IsInside(self, p, includeBoundary=True):
        if includeBoundary:
            if p.x >= self.lB.x and p.x <= self.tR.x:
                if p.y >= self.lB.y and p.y <= self.tR.y:
                    if p.z is None or (p.z >= self.lB.z and p.z <= self.tR.z):
                        # print("IsInside")
                        # sleep(0.05)
                        return True
        else:
            if p.x > self.lB.x and p.x < self.tR.x:
                if p.y > self.lB.y and p.y < self.tR.y:
                    if p.z is None or (p.z > self.lB.z and p.z < self.tR.z):
                        # print("IsInside")
                        # sleep(0.05)
                        return True
            # print("is Outside")
        return False
    def print(self):
        print("lB = (%lg, %lg, %lg)" % (self.lB.x, self.lB.y, self.lB.z))
        print("tR = (%lg, %lg, %lg)" % (self.tR.x, self.tR.y, self.tR.z))

class Sphere(Region):
    """docstring for Sphere."""
    mid = Point(0., 0.)
    radius = 0.

    def __init__(self, a, r):
        super(Region, self).__init__()
        self.mid = a
        self.radius = r

    def set_Geometry(self, a, r):
        self.mid = a
        self.radius = r
        pass

    def IsInside(self, p):
        if distance(p, self.mid) < self.radius:
            return True
        return False

    def print(self):
        print("mid = (%lg, %lg, %lg), radius = %lg" % (self.mid.x, self.mid.y, self.mid.z, self.radius))
