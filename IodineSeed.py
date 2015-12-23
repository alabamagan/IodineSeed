#!/usr/bin/python

import multiprocessing as mp
from scipy import integrate
import numpy as np

class IodineSeed(object):
    def __init__(self, center, length, orientation, activity):
        super(IodineSeed, self).__init__()
        self._center = np.array(center)
        self._length = length
        self._orientation = np.array(orientation)
        self._activity = activity

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
        if m_r == 0:
            return [0,0                       ]
        m_theta = m_theta % 360
        if m_theta > 180:
            m_theta = 360-m_theta
        m_directDisatance = m_r/np.sin(np.deg2rad(m_theta))
        m_edgeLeft = self._cosineLawGetLength(self._length/2, m_directDisatance, np.deg2rad(180. - m_theta))
        m_edgeRight = self._cosineLawGetLength(self._length/2, m_directDisatance, np.deg2rad(m_theta))
        m_degreeLeft = np.rad2deg(self._sinLawGetAngle(m_edgeLeft, np.deg2rad(180. - m_theta), m_directDisatance))
        m_degreeRight = 180. - np.rad2deg(self._sinLawGetAngle(m_edgeRight, np.deg2rad(m_theta), m_directDisatance))


        # print m_directDisatance, m_edgeLeft, m_edgeRight, m_degreeLeft, m_degreeRight
        m_intensityF = lambda theta: self._activity/np.sqrt(self._sinLawGetLength(m_directDisatance, np.deg2rad(m_degreeLeft), np.deg2rad(theta)))
        m_intensity = integrate.quad(m_intensityF, m_degreeLeft, m_degreeRight)
        print m_intensity
        return m_intensity

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

        m_scalarField = np.zeros([m_gridLengthX, m_gridLengthY, m_gridLengthZ], dtype=np.float32)
        for i in xrange(m_gridLengthX):
            for j in xrange(m_gridLengthY):
                for k in xrange(m_gridLengthZ):
                    l_l_l_coord = np.array([-self._boundX/2. + self._spacing[0]*i, -self._boundY/2. + self._spacing[1]*j, -self._boundZ/2. + self._spacing[2]*k])
                    l_l_l_relativeVect = l_l_l_coord - self._center
                    l_l_l_relativeUnitVect = l_l_l_relativeVect/np.linalg.norm(l_l_l_relativeVect)
                    l_l_l_theta = np.arcsin(np.dot(l_l_l_relativeUnitVect, self._orientation))
                    l_l_l_R = np.linalg.norm(l_l_l_relativeVect)*np.cos(l_l_l_theta)
                    m_scalarField[i,j,k] = self._energy(l_l_l_R, l_l_l_theta)[0]
                    print i,j,k

        # m_scalarField[m_gridLengthX:m_gridLengthX*2-1, m_gridLengthY:m_gridLengthY*2-1, ] =  m_scalarField[0:m_gridLengthX, 0:m_gridLengthY, 0:m_gridLengthZ]
        return m_scalarField

