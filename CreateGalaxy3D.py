from random import triangular, vonmisesvariate, choice
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import cos, sin, pi
import numpy as np

class Galaxy:
    # default_RANGE_Z = 1000.0 - limit of distribution of stars along the Z axis
    # default_R_GALAXY = 10000.0 - galaxy radius
    # default_KOL_STARS = 5000 - number of stars
    # default_START_EXPANSION = 30 - the beginning of the expansion of the galaxy
    # default_DELTA_X = 0 - axis offset (zero coordinate)
    # default_DELTA_Y = 0 - axis offset (zero coordinate)
    # default_DELTA_Z = 0 - axis offset (zero coordinate)
    RANGE_Z = 1000.0
    R_GALAXY = 10000.0
    KOL_STARS = 5000
    START_EXPANSION = 30
    DELTA_X = 0
    DELTA_Y = 0
    DELTA_Z = 0

    def __init__(self, name_galaxy = 'Milky_Way', radius_galaxy = R_GALAXY, range_z = RANGE_Z, start_expansion = START_EXPANSION, quantity_of_stars = KOL_STARS, delta_x = DELTA_X, delta_y = DELTA_Y, delta_z = DELTA_Z):
        self.name_galaxy = name_galaxy
        self.radius_galaxy = radius_galaxy
        self.range_z = range_z
        self.start_expansion = start_expansion
        self.quantity_of_stars = quantity_of_stars
        self.polsys_distances = []
        self.polsys_angles = []
        self.delta_x = delta_x
        self.delta_y = delta_y
        self.delta_z = delta_z
        for i in range(self.quantity_of_stars):
            self.polsys_distances.append(triangular(0, self.radius_galaxy))
            self.polsys_angles.append(vonmisesvariate(0, 0))
        self.X = self.transform_to_a_Cartesian_coordinate_system_XY('x')
        self.Y = self.transform_to_a_Cartesian_coordinate_system_XY('y')
        self.Z = self.point_sampling_Z()

    def transform_to_a_Cartesian_coordinate_system_XY(self, coordinate_str):
        # Convert polar coordinates to cartesian.
        # Преобразование полярных координат в декартовы.
        coordinaty_res = []
        if coordinate_str.lower() == 'x':
            for i in range(self.quantity_of_stars):
                coordinaty_res.append(self.polsys_distances[i] * cos(self.polsys_angles[i]) + self.delta_x)
        elif coordinate_str.lower() == 'y':
            for i in range(self.quantity_of_stars):
                coordinaty_res.append(self.polsys_distances[i] * sin(self.polsys_angles[i]) + self.delta_y)
        return coordinaty_res

    def point_sampling_Z(self):
        # Sample points on Z. Depending on the coefficient of the radius, a random or specified distribution by function of Z is specified. 
        # Выборка точек по Z. В зависимсости от коэффициента радиуса задается случайное или заданное функцией распределение по Z.
        z = []
        for j in range(self.quantity_of_stars):
            coeff_radius = (self.polsys_distances[j] * 100) / self.radius_galaxy
            if  coeff_radius >= self.start_expansion:
                z.append(triangular(-self.range_z, self.range_z) + self.delta_z)
            else:
                z.append(self.up_down_range_z(coeff_radius))
        return z

    def up_down_range_z(self, coeff_radius, coeff_expansion1 = 2.1, coeff_expansion2 = 1.9):
        # A special method for distributing stars around a supermassive black hole.
        # Специальный метод для распределения звезд вокруг сверхмассивной черной дыры.
        alpha = ((self.start_expansion - coeff_radius) / self.start_expansion) * pi/2
        z = (1 + (coeff_expansion1 * sin(alpha))) * triangular(coeff_expansion2 * sin(alpha), self.range_z) * choice([-1, 1])
        return z + self.delta_z

    def rotate_galaxy(self, around_axis, rotation_angle):
        # Rotate the galaxy at a given angle around a given axis.
        # Поворот галактики на заданный угол вокруг заданной оси.
        if around_axis.lower() == 'x':
            rotation_matrix = np.array([[1, 0, 0, 0], [0, cos(rotation_angle), sin(rotation_angle), 0], [0, -sin(rotation_angle), cos(rotation_angle), 0], [0, 0, 0, 1]])

        elif around_axis.lower() == 'y':
            rotation_matrix = np.array([[cos(rotation_angle), 0, -sin(rotation_angle), 0], [0, 1, 0, 0], [sin(rotation_angle), 0, cos(rotation_angle), 0], [0, 0, 0, 1]])
            
        for i in range(self.quantity_of_stars):
            position_star = np.array([self.X[i], self.Y[i], self.Z[i], 1])
            new_position_star = np.dot(position_star, rotation_matrix)
            self.X[i] = new_position_star[0]
            self.Y[i] = new_position_star[1]
            self.Z[i] = new_position_star[2]

    @classmethod
    def build_a_graph(cls, *objects):
        # Plotting galaxies
        # Построение графика галактик
        objects[0].X.append(30000)
        objects[0].Y.append(30000)
        objects[0].Z.append(30000)
        objects[0].X.append(-30000)
        objects[0].Y.append(-30000)
        objects[0].Z.append(-30000)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        for obj in objects:
            ax.scatter(obj.X, obj.Y, obj.Z, marker = '.')
        ax.set_xlabel('My X')
        ax.set_ylabel('My Y')
        ax.set_zlabel('My Z')
        plt.show()

# Create instances of our galaxies
# Создаем экземпляры наших галактик
milky_way = Galaxy('Milky_Way')
andromeda = Galaxy('Andromeda', radius_galaxy = 20000, start_expansion = 40, quantity_of_stars = 10000, delta_x = triangular(-9000, 9000), delta_y = triangular(-10000, 10000), delta_z = triangular(-5000, 5000))
andromeda.rotate_galaxy('x', pi / 6)

# We plot the collision of two galaxies
# Строим график столкновения двух галактик
Galaxy.build_a_graph(milky_way, andromeda)
