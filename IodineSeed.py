#!/usr/bin/python

import multiprocessing as mp
from scipy import integrate, interpolate
import numpy as np
import multiprocessing as mp

# ===========================================================
# Model 6711 I125 Amersham Health Model 6711
# ===========================================================
# 2D Anistropy Function
r_thetaTable = [[0.333, 0.370, 0.442, 0.488, 0.520, 0.550],
                [0.400, 0.429, 0.497, 0.535, 0.561, 0.587],
                [0.519, 0.537, 0.580, 0.609, 0.630, 0.645],
                [0.716, 0.705, 0.727, 0.743, 0.752, 0.760],
                [0.846, 0.834, 0.842, 0.846, 0.848, 0.852],
                [0.926, 0.925, 0.926, 0.926, 0.928, 0.928],
                [0.972, 0.972, 0.970, 0.969, 0.969, 0.969],
                [0.991, 0.991, 0.987, 0.987, 0.987, 0.987],
                [0.996, 0.996, 0.996, 0.995, 0.995, 0.995],
                [1.000, 1.000, 1.000, 0.999, 0.999, 0.999]]
theta = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80]
r = [0.5, 1, 2, 3, 4, 5]

# Line Source Dose-rate Constant cGy/(hU)
doseRateConstant = 0.965

# Line Source Radial-Dose Function
r_radialDose = [0.1, 0.15, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 5, 6, 7, 8, 9, 10]
radialDoseTable = [1.055, 1.078, 1.082, 1.071, 1.042, 1, 0.908, 0.814, 0.632, 0.496, 0.364, 0.270, 0.199, 0.148, 0.109,
                   0.0803]

anistropy = [0.973, 0.944, 0.941, 0.942, 0.943, 0.944]

# BackupPlan
doseTable = [[3, 0.965, 0.197, 0.0682, 0.0301, 0.0142, 0.00729, 0.00395],
             [1.524, 0.754, 0.183, 0.0657, 0.0294, 0.0139, 0.00720, 0.00391],
             [0.368, 0.427, 0.148, 0.0590, 0.0274, 0.0132, 0.00693, 0.00380],
             [0.161, 0.223, 0.110, 0.0499, 0.0244, 0.0122, 0.00651, 0.00361],
             [0.0879, 0.121, 0.0764, 0.0403, 0.0209, 0.0109, 0.00598, 0.00338],
             [0.0525, 0.0700, 0.0526, 0.0315, 0.0172, 0.00950, 0.00536, 0.00310],
             [0.0334, 0.0430, 0.0361, 0.0238, 0.0138, 0.00805, 0.00470, 0.00279],
             [0.0226, 0.0279, 0.0249, 0.0177, 0.0111, 0.00671, 0.00403, 0.00246],
             [0.0157, 0.0187, 0.0172, 0.0128, 0.00871, 0.00556, 0.00345, 0.00215],
             [0.0111, 0.0127, 0.0119, 0.00954, 0.00673, 0.00451, 0.00292, 0.00187],
             [0.00779, 0.00878, 0.00857, 0.00705, 0.00526, 0.00360, 0.00244, 0.00160],
             [0.00576, 0.00631, 0.00618, 0.00526, 0.00404, 0.00292, 0.00201, 0.00135],
             [0.00423, 0.00456, 0.00453, 0.00394, 0.00312, 0.00233, 0.00165, 0.00114],
             [0.00321, 0.00339, 0.00334, 0.00296, 0.00243, 0.00185, 0.00134, 0.000950],
             [0.00241, 0.00252, 0.00250, 0.00227, 0.00188, 0.00147, 0.00109, 0.000783]]
y = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0]
z = [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]


# ===========================================================
# End
# ===========================================================


class IodineSeed(object):
    def __init__(self, center, length, orientation, activity):
        super(IodineSeed, self).__init__()
        self._center = np.array(center)
        self._length = length
        self._orientation = np.array(orientation)
        self._activity = activity
        self._r_thetaTable = np.array(r_thetaTable)

        # Construct interpolated anistropy function
        self._anistropy = np.array(anistropy)
        self._theta = np.array(theta)
        self._r = np.array(r)
        self._anistropyF = self._constructInterpolationTable()

        # Construct interpolated raidialDoseFunction
        self._radialDoseTable = np.array(radialDoseTable)
        self._r_radialDose = np.array(r_radialDose)
        self._radialDoseTableF = interpolate.interp1d(self._r_radialDose, self._radialDoseTable, "linear",
                                                      bounds_error=False, fill_value="1")

        # Construct interpolation for  backup
        self._y = np.array(y)
        self._z = np.array(z)
        self._doseTable = np.array(doseTable)
        self._doseTableF = self._constructInterpolationDoseMap()

    def _sinLawGetAngle(self, m_r1, m_theta1, m_r2):
        m_theta2 = np.arcsin(m_r2 * np.sin(m_theta1) / m_r1)
        return m_theta2

    def _sinLawGetLength(self, m_r1, m_theta1, m_theta2):
        """
        Get Length using sine law

        :param m_r1:
        :param m_theta1:    [float] Angle 1 in radian
        :param m_theta2:    [float] Angle 2 in radian
        :return:
        """
        if m_theta1 == m_theta2:
            print "[Warning] Angle recieved equals."
            m_r2 = m_r1
        else:
            m_r2 = m_r1 * np.sin(m_theta2) / np.sin(m_theta1)
        return m_r2

    def _cosineLawGetLength(self, m_r1, m_r2, m_includeAngle):
        """
        Cosine law implementation

        :param m_r1:            [float]
        :param m_r2:            [float]
        :param m_includeAngle:  [float] Include Angle in radian
        :return:
        """
        m_out = np.sqrt(m_r1 ** 2 + m_r2 ** 2 - 2 * m_r1 * m_r2 * np.cos(m_includeAngle))
        return m_out

    def _energy(self, m_r, m_theta):
        m_r = np.abs(m_r)
        if m_r < 0.1:
            return [0, 0]
            # m_theta = m_theta % 360
            # if m_theta > 180:
            # m_theta = 360-m_theta
        m_directDisatance = m_r / np.sin(np.deg2rad(m_theta))
        m_edgeLeft = self._cosineLawGetLength(self._length / 2, m_directDisatance, np.deg2rad(180. - m_theta))
        m_edgeRight = self._cosineLawGetLength(self._length / 2, m_directDisatance, np.deg2rad(m_theta))
        m_degreeLeft = np.rad2deg(self._sinLawGetAngle(m_edgeLeft, np.deg2rad(180. - m_theta), m_directDisatance))
        m_degreeRight = np.rad2deg(self._sinLawGetAngle(m_edgeRight, np.deg2rad(m_theta), m_directDisatance))
        print m_degreeLeft, m_degreeRight

        # print m_directDisatance, m_edgeLeft, m_edgeRight, m_degreeLeft, m_degreeRight
        m_intensityF = lambda base: self._activity / self._cosineLawGetLength(m_edgeLeft, base,
                                                                              np.deg2rad(m_degreeLeft)) ** 2
        m_intensity = integrate.quad(m_intensityF, 0, self._length)
        print m_intensity
        return m_intensity

    def _lineSourceGL(self, m_r, m_theta):
        m_edgeLeft = self._cosineLawGetLength(self._length / 2, m_r, np.deg2rad(180. - m_theta))
        m_edgeRight = self._cosineLawGetLength(self._length / 2, m_r, np.deg2rad(m_theta))
        m_degreeLeft = np.rad2deg(self._sinLawGetAngle(m_edgeLeft, np.deg2rad(180. - m_theta), m_r))
        m_degreeRight = np.rad2deg(self._sinLawGetAngle(m_edgeRight, np.deg2rad(m_theta), m_r))
        m_degreeBeta = np.abs(m_degreeLeft - m_degreeRight)
        print m_degreeLeft, m_degreeRight, m_degreeBeta, m_edgeLeft, m_edgeRight
        if m_theta != 0:
            m_val = np.deg2rad(m_degreeBeta) / (self._length * m_r * np.sin(np.deg2rad(m_theta)))
        else:
            m_val = 1 / np.abs(m_r ** 2 - self._length ** 2 / 4.)
        return m_val

    def _constructInterpolationTable(self):
        m_meshgrid = np.meshgrid(self._r, self._theta)
        m_r = m_meshgrid[0].flatten()
        m_theta = m_meshgrid[1].flatten()
        m_table = self._r_thetaTable.flatten()
        m_f = interpolate.interp2d(m_r, m_theta, m_table, kind="quintic")
        return m_f

    def _constructInterpolationDoseMap(self):
        m_meshgrid = np.meshgrid(self._y, self._z)
        m_y = m_meshgrid[0].flatten()
        m_z = m_meshgrid[1].flatten()
        m_table = self._doseTable.flatten()
        m_f = interpolate.interp2d(m_y, m_z, m_table, kind="quintic")
        return m_f

    def SetSize(self, x, y, z):
        self._boundX = x
        self._boundY = y
        self._boundZ = z
        pass

    def SetSpacing(self, x, y, z):
        self._spacing = [x, y, z]
        pass

    def Update(self):
        m_gridLengthX = np.int(self._boundX / self._spacing[0])
        m_gridLengthY = np.int(self._boundY / self._spacing[1])
        m_gridLengthZ = np.int(self._boundZ / self._spacing[2])

        m_CPUnumber = mp.cpu_count()
        pool = mp.Pool(m_CPUnumber)
        result = []
        m_scalarField = np.zeros([m_gridLengthX, m_gridLengthY, m_gridLengthZ], dtype=np.float32)
        for i in xrange(m_gridLengthX):
            for j in xrange(m_gridLengthY):
                for k in xrange(m_gridLengthZ):
                    l_l_l_coord = np.array(
                        [-self._boundX / 2. + self._spacing[0] * i, -self._boundY / 2. + self._spacing[1] * j,
                         -self._boundZ / 2. + self._spacing[2] * k])
                    l_l_l_relativeVect = l_l_l_coord - self._center
                    l_l_l_relativeUnitVect = l_l_l_relativeVect / np.linalg.norm(l_l_l_relativeVect)
                    l_l_l_theta = np.rad2deg(np.arccos(np.dot(l_l_l_relativeUnitVect, self._orientation)))
                    l_l_l_R = np.abs(np.linalg.norm(l_l_l_relativeVect) * np.sin(np.deg2rad(l_l_l_theta)))
                    # result.append(pool.apply_async(self._energy(l_l_l_R, l_l_l_theta)[0]))
                    m_scalarField[i, j, k] = self._mg(l_l_l_R, l_l_l_theta)
                    print i, j, k
        # pool.close()
        # pool.join()
        # print result
        # m_scalarField[m_gridLengthX:m_gridLengthX*2-1, m_gridLengthY:m_gridLengthY*2-1, ] =  m_scalarField[0:m_gridLengthX, 0:m_gridLengthY, 0:m_gridLengthZ]
        return m_scalarField

    def Update2D(self):
        m_gridLengthX = np.int(self._boundX / (2 * self._spacing[0]))
        m_gridLengthY = np.int(self._boundY / (2 * self._spacing[1]))

        print m_gridLengthX, m_gridLengthY

        m_scalarField = np.zeros([m_gridLengthX, m_gridLengthY], dtype=np.float32)
        for i in xrange(m_gridLengthX):
            for j in xrange(m_gridLengthY):
                l_l_l_coord = np.array([self._spacing[0] * i, self._spacing[1] * j])
                l_l_l_relativeVect = self._center[0:2] - l_l_l_coord
                l_l_l_relativeUnitVect = l_l_l_relativeVect / np.linalg.norm(l_l_l_relativeVect)
                l_l_l_R = np.linalg.norm(l_l_l_relativeVect)
                l_l_l_theta = np.rad2deg(np.arctan(l_l_l_coord[1]/l_l_l_coord[0]))
                m_scalarField[i, j] =  self._lineSourceGL(l_l_l_R, l_l_l_theta)
        return np.flipud(m_scalarField)

    def GenerateFromInterpolation2D(self):
        m_gridLengthX = np.int(self._boundX / self._spacing[0])
        m_gridLengthY = np.int(self._boundY / self._spacing[1])

        m_scalarField = np.zeros([m_gridLengthX, m_gridLengthY], dtype=np.float32)
        for i in xrange(m_gridLengthX):
            for j in xrange(m_gridLengthY):
                l_l_l_coord = np.array([self._spacing[0] * i, self._spacing[1] * j])
                l_l_l_relativeVect = l_l_l_coord - self._center[0:2]
                l_l_l_relativeUnitVect = l_l_l_relativeVect / np.linalg.norm(l_l_l_relativeVect)
                l_l_l_theta = np.rad2deg(np.arccos(np.dot(l_l_l_relativeUnitVect, self._orientation[0:2])))
                l_l_l_R = np.linalg.norm(l_l_l_relativeVect)
                # result.append(pool.apply_async(self._energy(l_l_l_R, l_l_l_theta)[0]))
                # m_scalarField[i,j] = self._lineSourceGL(l_l_l_R, l_l_l_theta)/self._lineSourceGL(1, 0) * self._radialDoseTableF(l_l_l_R) * self._anistropyF(l_l_l_R, l_l_l_theta)
                m_scalarField[i, j] = self._doseTableF(l_l_l_coord[0], l_l_l_coord[1])
                # print i,j
        return np.flipud(m_scalarField)
