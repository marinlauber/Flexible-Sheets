# -*- coding: utf-8 -*-

# Imports
import numpy as np
import vtk


def read_vti(fname):
	reader = vtk.vtkXMLPImageDataReader()
	reader.SetFileName(fname)
	reader.Update()
	data = reader.GetOutput()
	pointData = data.GetPointData()

	sh = data.GetDimensions()[::-1]
	ndims = len(sh)

	# get vector field
	v = np.array(pointData.GetVectors("Velocity")).reshape(sh + (ndims,))
	vec = []
	for d in range(ndims):
		a = v[..., d]
		vec.append(a)
	# get scalar field
	sca = np.array(pointData.GetScalars('Pressure')).reshape(sh + (1,))

	# Generate grid
	# nPoints = data.GetNumberOfPoints()
	(xmin, xmax, ymin, ymax, zmin, zmax) = data.GetBounds()
	grid3D = np.mgrid[xmin:xmax + 1, ymin:ymax + 1, zmin:zmax + 1]

	return np.transpose(np.array(vec), (0,3,2,1)), np.transpose(sca, (0,3,2,1)), grid3D


def read_vtr(fname):
	reader = vtk.vtkXMLPRectilinearGridReader()
	reader.SetFileName(fname)
	reader.Update()
	data = reader.GetOutput()
	pointData = data.GetPointData()

	sh = data.GetDimensions()[::-1]
	ndims = len(sh)

	# get vector field
	v = np.array(pointData.GetVectors("Velocity")).reshape(sh + (ndims,))
	vec = []
	for d in range(ndims):
		a = v[..., d]
		vec.append(a)
	vec = np.array(vec)

	# get scalar field
	sca = np.array(pointData.GetScalars('Pressure')).reshape(sh + (1,))

	# get grid
	x = np.array(data.GetXCoordinates())
	y = np.array(data.GetYCoordinates())
	z = np.array(data.GetZCoordinates())

	return np.transpose(vec, (0,3,2,1)), np.transpose(sca, (3,2,1,0))[0,:,:,:], np.array([x, y, z], dtype=object)
