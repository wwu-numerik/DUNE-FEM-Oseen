#!/usr/env/python
# -*- coding: utf-8 -*-

class DimensionIncompatibleException(Exception):
	pass

class ImpossibleException(Exception):
	pass

import math
from euclid import *

class PLCPointList:
	"""static var numbering all vertices"""
	global_vertices = []

	def __init__(self, dim ):
		self.dim = dim
		self.verts = dict()

	def appendVert(self,x):
		if not x in PLCPointList.global_vertices:
			PLCPointList.global_vertices.append(x)
			glob_idx = len(PLCPointList.global_vertices)
		else:
			glob_idx = PLCPointList.global_vertices.index(x)
		self.verts[glob_idx] = x
		return glob_idx

	def simpleString(self):
		ret = ''
		for idx,p in self.verts.iteritems():
			ret += 'Vertix %4d \t%s\n'% (idx,str(p))
		return ret

class Simplex:
	def __init__(self,v1,v2,v3):
		self.v1 = v1
		self.v2 = v2
		self.v3 = v3

	def __repr__(self):
		return 'Simplex (%4d,%4d,%4d)'%(self.v1,self.v2,self.v3)

class FanningSimplexList:
	def __init__(self,center_idx,bid):
		self.center_idx = center_idx
		self.simplices 	= []
		self.boundaryId = bid
		self.vertex_idx = []

	def addVertex(self,v_idx):
		self.vertex_idx.append( v_idx )
		if len(self.vertex_idx) > 1:
			self.simplices.append( Simplex( self.center_idx, self.vertex_idx[-1], self.vertex_idx[-2]  ) )

	def close(self):
		"""finish last simplex"""
		self.simplices.append( Simplex( self.center_idx, self.vertex_idx[-1], self.vertex_idx[0]  ) )

	def __repr__(self):
		ret = 'FanningSimplexList  for boundary ID %d\n'%(self.boundaryId)
		i = 0
		for s in self.simplices:
			ret += '%4d %s\n'%(i,s)
			i += 1
		return ret

	def __str__(self):
		return self.__repr__()

class FullGrid:
	def __init__(self,f1,default_Bid):
		self.fans 					= [f1]
		self.default_Bid 			= default_Bid
		self.connecting_simplices 	= []

	def connect(self,new_f):
		f1 = self.fans[-1]
		f2 = new_f
		if len(f1.vertex_idx) !=  len(f2.vertex_idx):
			raise DimensionIncompatibleException()
		b_len = len(f1.vertex_idx)
		for i in range ( 0, b_len  ):
			self.connecting_simplices.append( \
				Simplex(	f1.vertex_idx[i-1],
							f2.vertex_idx[i-1], \
							f1.vertex_idx[i] ) )
			self.connecting_simplices.append( \
				Simplex(	f1.vertex_idx[i],
							f2.vertex_idx[i], \
							f2.vertex_idx[i-1] ) )
	def __str__(self):
		ret = 'connecting simplices %d\n'%(len(self.connecting_simplices))
		for s in self.connecting_simplices:
			ret += str(s) + '\n'
		return ret

	def outputPLC(self,fn):
		out = None
		try:
			out = open(fn,'w')
		except:
			raise ImpossibleException()
		out.write( '%d 3 0 %d\n'%(len(PLCPointList.global_vertices),3) )#3 bids
		cVert = 1
		for v in PLCPointList.global_vertices:
				out.write( '%d %f %f %f\n'%(cVert,v.x,v.y,v.z) )
				cVert += 1
		boundary_segment_count = len(self.connecting_simplices)
		for f in self.fans:
			boundary_segment_count += len(f.simplices)
		out.write( '%d 1\n'%(boundary_segment_count) )
		cBseg = 0
		for f in self.fans:
			for b1 in f.simplices:
				out.write( '%d %d %d %d %d\n'%(3,b1.v1,b1.v2,b1.v3,f.boundaryId)  )
				cBseg += 1
		for c in self.connecting_simplices:
			out.write( '%d %d %d %d %d\n'%(3,c.v1,c.v2,c.v3,self.default_Bid)  )
			cBseg += 1
		out.write( '%d\n'%(0))










