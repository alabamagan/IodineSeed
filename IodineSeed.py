#!/usr/bin/python

import multiprocessing as mp
from scipy import integrate, interpolate
import numpy as np
import multiprocessing as mp

r_thetaTable = [[0.333,0.370,0.442,0.488,0.520,0.550],
                [0.400,0.429,0.497,0.535,0.561,0.587],
                [0.519,0.537,0.580,0.609,0.630,0.645],
                [0.716,0.705,0.727,0.743,0.752,0.760],
                [0.846,0.834,0.842,0.846,0.848,0.852],
                [0.926,0.925,0.926,0.926,0.928,0.928],
                [0.972,0.972,0.970,0.969,0.969,0.969],
                [0.991,0.991,0.987,0.987,0.987,0.987],
                [0.996,0.996,0.996,0.995,0.995,0.995],
                [1.000,1.000,1.000,0.999,0.999,0.999]]

theta = [0,5,10,20,30,40,50,60,70,80]
r = [0.5,1,2,3,4,5]

anistropy = [0.973,0.944,0.941,0.942,0.943,0.944]

class IodineSeed(object):
    def __init__(self, center, length, orientation, activity):
        super(IodineSeed, self).__init__()
        self._center = np.array(center)
        self._length = length
        self._orientation = np.array(orientation)
        self._activity = activity
        self._r_thetaTable = np.array(r_thetaTable)
        self._anistropy = np.array(anistropy)
        self._theta = np.array(theta)
        self._r = np.array(r)
        self._anistropyF = self._constructInterpolationTable()


    def _sinLawGetAngle(self, m_r1, m_theta1, m_r2):
        m_theta2 = np.arcsin(m_r2*np.sin(m_theta1)/m_r1)
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
            m_r2 = m_r1*np.sin(m_theta2)/np.sin(m_theta1)
        return m_r2

    def _cosineLawGetLength(self, m_r1, m_r2, m_includeAngle):
        """
        Cosine law implementation

        :param m_r1:            [float]
        :param m_r2:            [float]
        :param m_includeAngle:  [float] Include Angle in radian
        :return:
        """
        m_out = np.sqrt(m_r1**2 + m_r2**2 - 2*m_r1*m_r2*np.cos(m_includeAngle))
        return m_out

    def _energy(self, m_r, m_theta):
        m_r = np.abs(m_r)
        if m_r < 0.1:
            return [0,0]
        # m_theta = m_theta % 360
        # if m_theta > 180:
            # m_theta = 360-m_theta
        m_directDisatance = m_r/np.sin(np.deg2rad(m_theta))
        m_edgeLeft = self._cosineLawGetLength(self._length/2, m_directDisatance, np.deg2rad(180. - m_theta))
        m_edgeRight = self._cosineLawGetLength(self._length/2, m_directDisatance, np.deg2rad(m_theta))
        m_degreeLeft = np.rad2deg(self._sinLawGetAngle(m_edgeLeft, np.deg2rad(180. - m_theta), m_directDisatance))
        m_degreeRight = np.rad2deg(self._sinLawGetAngle(m_edgeRight, np.deg2rad(m_theta), m_directDisatance))
        print m_degreeLeft, m_degreeRight

        # print m_directDisatance, m_edgeLeft, m_edgeRight, m_degreeLeft, m_degreeRight
        m_intensityF = lambda base: self._activity/self._cosineLawGetLength(m_edgeLeft,base,np.deg2rad(m_degreeLeft))**2
        m_intensity = integrate.quad(m_intensityF, 0, self._length)
        print m_intensity
        return m_intensity


    def _lineSourceGL(self,m_r, m_theta):

        if m_theta == 0:


    def _constructInterpolationTable(self):
        m_meshgrid = np.meshgrid(self._r,self._theta)
        m_r = m_meshgrid[0].flatten()
        m_theta = m_meshgrid[1].flatten()
        m_table = self._r_thetaTable.flatten()
        m_f = interpolate.interp2d(m_r,m_theta,m_table, kind="quintic")
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
        m_gridLengthX = np.int(self._boundX/self._spacing[0])
        m_gridLengthY = np.int(self._boundY/self._spacing[1])
        m_gridLengthZ = np.int(self._boundZ/self._spacing[2])

        m_CPUnumber = mp.cpu_count()
        pool = mp.Pool(m_CPUnumber)
        result = []
        m_scalarField = np.zeros([m_gridLengthX, m_gridLengthY, m_gridLengthZ], dtype=np.float32)
        for i in xrange(m_gridLengthX):
            for j in xrange(m_gridLengthY):
                for k in xrange(m_gridLengthZ):
                    l_l_l_coord = np.array([-self._boundX/2. + self._spacing[0]*i, -self._boundY/2. + self._spacing[1]*j, -self._boundZ/2. + self._spacing[2]*k])
                    l_l_l_relativeVect = l_l_l_coord - self._center
                    l_l_l_relativeUnitVect = l_l_l_relativeVect/np.linalg.norm(l_l_l_relativeVect)
                    l_l_l_theta = np.rad2deg(np.arccos(np.dot(l_l_l_relativeUnitVect, self._orientation)))
                    l_l_l_R = np.abs(np.linalg.norm(l_l_l_relativeVect)*np.sin(np.deg2rad(l_l_l_theta)))
                    # result.append(pool.apply_async(self._energy(l_l_l_R, l_l_l_theta)[0]))
                    m_scalarField[i,j,k] = self._mg(l_l_l_R, l_l_l_theta)
                    print i,j,k
        # pool.close()
        # pool.join()
        # print result
        # m_scalarField[m_gridLengthX:m_gridLengthX*2-1, m_gridLengthY:m_gridLengthY*2-1, ] =  m_scalarField[0:m_gridLengthX, 0:m_gridLengthY, 0:m_gridLengthZ]
        return m_scalarField

    def Update2D(self):
        m_gridLengthX = np.int(self._boundX/self._spacing[0])
        m_gridLengthY = np.int(self._boundY/self._spacing[1])

        m_CPUnumber = mp.cpu_count()
        pool = mp.Pool(m_CPUnumber)
        result = []
        m_scalarField = np.zeros([m_gridLengthX/2, m_gridLengthY/2], dtype=np.float32)
        for i in xrange(m_gridLengthX/2):
            for j in xrange(m_gridLengthY/2):
                l_l_l_coord = np.array([self._spacing[0]*i, self._spacing[1]*j])
                l_l_l_relativeVect = l_l_l_coord - self._center[0:2]
                l_l_l_relativeUnitVect = l_l_l_relativeVect/np.linalg.norm(l_l_l_relativeVect)
                l_l_l_theta = np.rad2deg(np.arccos(np.dot(l_l_l_relativeUnitVect, self._orientation[0:2])))
                l_l_l_R = np.linalg.norm(l_l_l_relativeVect)
                # result.append(pool.apply_async(self._energy(l_l_l_R, l_l_l_theta)[0]))
                m_scalarField[i,j] = self._anistropyF(l_l_l_R, l_l_l_theta)
                print i,j
        # pool.close()
        # pool.join()
        # print result
        # m_scalarField[m_gridLengthX:m_gridLengthX*2-1, m_gridLengthY:m_gridLengthY*2-1, ] =  m_scalarField[0:m_gridLengthX, 0:m_gridLengthY, 0:m_gridLengthZ]
        return np.flipud(m_scalarField)

