
''' The main window of the curver GUI application. '''

import curver
import curver.application

import re
import os
import io
import sys
import pickle
from math import sin, cos, pi, ceil, sqrt
from itertools import combinations

try:
	import Tkinter as TK
	import tkFileDialog
	import tkMessageBox
except ImportError:  # Python 3.
	try:
		import tkinter as TK
		import tkinter.filedialog as tkFileDialog
		import tkinter.messagebox as tkMessageBox
	except ImportError:
		raise ImportError('Tkinter not available.')

try:
	import ttk as TTK
except ImportError:  # Python 3.
	try:
		from tkinter import ttk as TTK
	except ImportError:
		raise ImportError('Ttk not available.')

# Some constants.
if sys.platform in ['darwin']:
	COMMAND = {
		'close': 'Command+W',
		}
	COMMAND_KEY = {
		'close': '<Command-w>',
		}
else:
	COMMAND = {
		'close': 'Ctrl+W',
		}
	COMMAND_KEY = {
		'close': '<Control-w>',
		}

# Vectors to offset a label by to produce backing.
OFFSETS = [(1.5*cos(2 * pi * i / 12), 1.5*sin(2 * pi * i / 12)) for i in range(12)]

# Colours of things.
DEFAULT_EDGE_LABEL_COLOUR = 'black'
DEFAULT_EDGE_LABEL_BG_COLOUR = 'white'
MAX_DRAWABLE = 1000  # Maximum weight of a multicurve to draw fully.

def dot(a, b):
	return a[0] * b[0] + a[1] * b[1]

class CurverApplication(object):
	def __init__(self, parent):
		self.parent = parent
		self.options = curver.application.Options(self)
		self.colour_picker = curver.application.ColourPalette()
		
		###
		self.canvas = TK.Canvas(self.parent, height=1, bg='#dcecff', takefocus=True)
		self.canvas.pack(padx=6, pady=6, fill='both', expand=True)
		self.canvas.bind('<Button-1>', self.canvas_left_click)
		self.canvas.bind('<Button-3>', self.canvas_right_click)
		self.canvas.bind('<Motion>', self.canvas_move)
		self.canvas.bind('<FocusOut>', self.canvas_focus_lost)
		###
		
		# Create the menus.
		# Make sure to start the Lamination and Mapping class menus disabled.
		self.menubar = TK.Menu(self.parent)
		app_font = self.options.application_font  # Get a shorter name.
		
		self.filemenu = TK.Menu(self.menubar, tearoff=0)
		self.filemenu.add_command(label='Export image...', command=self.export_image)
		self.filemenu.add_separator()
		self.filemenu.add_command(label='Exit', command=self.quit, accelerator=COMMAND['close'])
		self.menubar.add_cascade(label='File', menu=self.filemenu)
		
		##########################################
		self.settingsmenu = TK.Menu(self.menubar, tearoff=0)
		
		self.sizemenu = TK.Menu(self.menubar, tearoff=0)
		self.sizemenu.add_radiobutton(label='Small', var=self.options.size_var, value=curver.application.options.SIZE_SMALL)
		self.sizemenu.add_radiobutton(label='Medium', var=self.options.size_var, value=curver.application.options.SIZE_MEDIUM)
		self.sizemenu.add_radiobutton(label='Large', var=self.options.size_var, value=curver.application.options.SIZE_LARGE)
		# self.sizemenu.add_radiobutton(label='Extra large', var=self.options.size_var, value=curver.application.options.SIZE_XLARGE)
		
		self.edgelabelmenu = TK.Menu(self.menubar, tearoff=0)
		self.edgelabelmenu.add_radiobutton(label=curver.application.options.LABEL_EDGES_NONE, var=self.options.label_edges_var)
		self.edgelabelmenu.add_radiobutton(label=curver.application.options.LABEL_EDGES_INDEX, var=self.options.label_edges_var)
		self.edgelabelmenu.add_radiobutton(label=curver.application.options.LABEL_EDGES_GEOMETRIC, var=self.options.label_edges_var)
		self.edgelabelmenu.add_radiobutton(label=curver.application.options.LABEL_EDGES_ALGEBRAIC, var=self.options.label_edges_var)
		self.edgelabelmenu.add_separator()
		self.edgelabelmenu.add_checkbutton(label='Projectivise', var=self.options.projectivise_var)
		
		self.zoommenu = TK.Menu(self.menubar, tearoff=0)
		self.zoommenu.add_command(label='Zoom in', command=self.zoom_in, accelerator='+')
		self.zoommenu.add_command(label='Zoom out', command=self.zoom_out, accelerator='-')
		self.zoommenu.add_command(label='Zoom to drawing', command=self.zoom_to_drawing, accelerator='0')
		
		self.settingsmenu.add_cascade(label='Sizes', menu=self.sizemenu)
		self.settingsmenu.add_cascade(label='Edge label', menu=self.edgelabelmenu)
		self.settingsmenu.add_cascade(label='Zoom', menu=self.zoommenu)
		self.settingsmenu.add_checkbutton(label='Show internal edges', var=self.options.show_internals_var)
		self.settingsmenu.add_checkbutton(label='Show edge orientations', var=self.options.show_orientations_var)
		
		self.menubar.add_cascade(label='Settings', menu=self.settingsmenu)
		
		self.helpmenu = TK.Menu(self.menubar, tearoff=0)
		self.helpmenu.add_command(label='Help', command=self.show_help, accelerator='F1')
		self.helpmenu.add_separator()
		self.helpmenu.add_command(label='About', command=self.show_about)
		
		self.menubar.add_cascade(label='Help', menu=self.helpmenu)
		self.parent.config(menu=self.menubar)
		
		self.parent.bind(COMMAND_KEY['close'], lambda event: self.quit())
		self.parent.bind('<Key>', self.parent_key_press)
		
		self.parent.protocol('WM_DELETE_WINDOW', self.quit)
		
		self.vertices = []
		self.edges = []
		self.triangles = []
		self.curve_components = []
		
		self.zeta = 0
		
		self.output = None
	
	def initialise(self):
		self.destroy_all_vertices()
		self.colour_picker.reset()
		
		self.unsaved_work = False
		return True
	
	def load(self, load_from=None):
		''' Load up some information.
		
		We can load from:
		'''
		
		try:
			if isinstance(load_from, curver.kernel.Triangulation):
				triangulation = load_from
				lamination = triangulation.empty_lamination()
			elif isinstance(load_from, curver.kernel.Lamination):
				triangulation = load_from.triangulation
				lamination = load_from
			
			# We don't have any vertices or edges, so we'll create a triangulation ourselves.
			vertices, edges = [], []
			
			# Get a dual tree.
			dual_tree = triangulation.dual_tree()
			components = triangulation.components()
			num_components = len(components)
			# Make sure we get the right sizes of things.
			self.parent.update_idletasks()
			w = int(self.canvas.winfo_width())
			h = int(self.canvas.winfo_height())
			
			# We will layout the components in a p x q grid.
			# Aim to maximise r === min(w / p, h / q) subject to pq >= num_components.
			# Note that there is probably a closed formula for the optimal value of p (and so q).
			p = max(range(1, num_components+1), key=lambda p: min(w / p, h / ceil(float(num_components) / p)))
			q = int(ceil(float(num_components) / p))
			
			r = min(w / p, h / q) * (1 + self.options.zoom_fraction) / 4
			dx = w / p
			dy = h / q
			
			num_used_vertices = 0
			for index, component in enumerate(components):
				# Get the number of triangles.
				n = len(component) // 3  # Remember component double counts edges.
				ngon = n + 2
				
				# Create the vertices.
				for i in range(ngon):
					vertices.append((
						dx * (index % p) + dx / 2 + r * sin(2 * pi * (i + 0.5) / ngon),
						dy * int(index / p) + dy / 2 + r * cos(2 * pi * (i + 0.5) / ngon)
						))
				
				def num_descendants(edge_label):
					''' Return the number of triangles that can be reached in the dual tree starting at the given edge_label. '''
					
					corner = triangulation.corner_of_edge(edge_label)
					left = (1 + sum(num_descendants(~(corner.labels[2])))) if dual_tree[corner.indices[2]] else 0
					right = (1 + sum(num_descendants(~(corner.labels[1])))) if dual_tree[corner.indices[1]] else 0
					
					return left, right
				
				initial_edge_index = min(i for i in component if i >= 0 and not dual_tree[i])
				to_extend = [(num_used_vertices, num_used_vertices+1, initial_edge_index)]
				# Hmmm, need to be more careful here to ensure that we correctly orient the edges.
				edges.append((num_used_vertices+1, num_used_vertices+0, initial_edge_index, None))
				while to_extend:
					source_vertex, target_vertex, label = to_extend.pop()
					left, right = num_descendants(label)
					far_vertex = target_vertex + left + 1
					corner = triangulation.corner_of_edge(label)
					
					if corner.labels[2] == corner.indices[2]:
						edges.append((far_vertex, target_vertex, corner.indices[2], None))
					else:
						edges.append((target_vertex, far_vertex, corner.indices[2], None))
					if corner.labels[1] == corner.indices[1]:
						edges.append((source_vertex, far_vertex, corner.indices[1], None))
					else:
						edges.append((far_vertex, source_vertex, corner.indices[1], None))
					
					if left > 0:
						to_extend.append((far_vertex, target_vertex, ~(corner.labels[2])))
					
					if right > 0:
						to_extend.append((source_vertex, far_vertex, ~(corner.labels[1])))
				num_used_vertices = len(vertices)
			
			# Glue together sides with the same index.
			for i, j in combinations(range(len(edges)), r=2):
				if edges[i][2] == edges[j][2]:
					edges[i] = (edges[i][0], edges[i][1], edges[i][2], j)
					edges[j] = (edges[j][0], edges[j][1], edges[j][2], i)
			
			if not self.initialise():
				return
			
			# Create the vertices.
			for vertex in vertices:
				self.create_vertex(vertex)
			
			# Create the edges.
			for edge in edges:
				start_index, end_index, edge_index, glued_to_index = edge
				self.create_edge(self.vertices[start_index], self.vertices[end_index])
			
			# Create the edge identifications.
			for index, edge in enumerate(edges):
				start_index, end_index, edge_index, glued_to_index = edge
				if glued_to_index is not None and glued_to_index > index:
					self.create_edge_identification(self.edges[index], self.edges[glued_to_index])
			
			# Set the correct edge indices.
			for index, edge in enumerate(edges):
				start_index, end_index, edge_index, glued_to_index = edge
				self.edges[index].index = edge_index
			
			self.equipped_triangulation = equipped_triangulation
			
			self.zoom_to_drawing()
			
			# Get the correct empty lamination.
			self.lamination_to_canvas(lamination)
		except (curver.AssumptionError, IndexError, ValueError) as error:
			tkMessageBox.showerror('Load Error', error.message)
	
	def export_image(self):
		path = tkFileDialog.asksaveasfilename(defaultextension='.ps', filetypes=[('postscript files', '.ps'), ('all files', '.*')], title='Export Image')
		if path:
			try:
				self.canvas.postscript(file=path, colormode='color')
			except IOError:
				tkMessageBox.showwarning('Export Error', 'Could not open: %s' % path)
	
	def quit(self):
		# Write down our current state for output. If we are incomplete then this is just None.
		
		if self.initialise():
			# Apparantly there are some problems with comboboxes, see:
			#  http://stackoverflow.com/questions/15448914/python-tkinter-ttk-combobox-throws-exception-on-quit
			self.parent.eval('::ttk::CancelRepeat')
			self.parent.destroy()
			self.parent.quit()
	
	def show_help(self):
		curver.doc.open_documentation()
	
	def show_about(self):
		tkMessageBox.showinfo('About', 'curver (Version %s).\nCopyright (c) Mark Bell 2013.' % curver.__version__)
	
	def translate(self, dx, dy):
		for vertex in self.vertices:
			vertex[0] = vertex[0] + dx
			vertex[1] = vertex[1] + dy
		
		for curve_component in self.curve_components + self.train_track_blocks:
			for i in range(len(curve_component.vertices)):
				curve_component.vertices[i] = curve_component.vertices[i][0] + dx, curve_component.vertices[i][1] + dy
		
		self.canvas.move('all', dx, dy)
	
	def zoom(self, scale):
		for vertex in self.vertices:
			vertex[0], vertex[1] = scale * vertex[0], scale * vertex[1]
			vertex.update()
		for edge in self.edges:
			edge.update()
		for triangle in self.triangles:
			triangle.update()
		for curve_component in self.curve_components + self.train_track_blocks:
			for i in range(len(curve_component.vertices)):
				curve_component.vertices[i] = scale * curve_component.vertices[i][0], scale * curve_component.vertices[i][1]
			curve_component.update()
		self.build_edge_labels()
		self.redraw()
	
	def zoom_in(self):
		self.zoom_centre(1.05)
	
	def zoom_out(self):
		self.zoom_centre(0.95)
	
	def zoom_centre(self, scale):
		self.parent.update_idletasks()
		cw = int(self.canvas.winfo_width())
		ch = int(self.canvas.winfo_height())
		self.translate(-cw / 2, -ch / 2)
		self.zoom(scale)
		self.translate(cw / 2, ch / 2)
	
	def zoom_to_drawing(self):
		self.parent.update_idletasks()
		box = self.canvas.bbox('all')
		if box is not None:
			x0, y0, x1, y1 = box
			cw = int(self.canvas.winfo_width())
			ch = int(self.canvas.winfo_height())
			cr = min(cw, ch)
			
			w, h = x1 - x0, y1 - y0
			r = max(w, h)
			
			self.translate(-x0 - w / 2, -y0 - h / 2)
			self.zoom(self.options.zoom_fraction * float(cr) / r)
			self.translate(cw / 2, ch / 2)
	
	
	######################################################################
	
	
	def create_vertex(self, point):
		self.vertices.append(curver.application.CanvasVertex(self.canvas, point, self.options))
		self.unsaved_work = True
		self.redraw()
		self.build_equipped_triangulation()
		return self.vertices[-1]
	
	def destroy_vertex(self, vertex=None):
		if vertex is None: vertex = self.vertices[-1]
		if self.selected_object == vertex: self.select_object(None)
		while True:
			for edge in self.edges:
				if edge[0] == vertex or edge[1] == vertex:
					self.destroy_edge(edge)
					break
			else:
				break
		self.canvas.delete(vertex.drawn)
		self.vertices.remove(vertex)
		self.unsaved_work = True
		self.redraw()
		self.build_equipped_triangulation()
	
	def destroy_all_vertices(self):
		while self.vertices:
			self.destroy_vertex()
		self.build_equipped_triangulation()
	
	def create_edge(self, v1, v2):
		# Check that the vertices are distinct
		if len(set([v1, v2])) != 2:
			return None
		
		# Check that this edge doesn't already exist.
		if any(set([edge[0], edge[1]]) == set([v1, v2]) for edge in self.edges):
			return None
		
		# Check that this edge doesn't intersect an existing one.
		if any(curver.application.lines_intersect(edge[0], edge[1], v1, v2, self.options.float_error, True)[1] for edge in self.edges):
			return None
		
		e0 = curver.application.CanvasEdge(self.canvas, [v1, v2], self.options)
		self.edges.append(e0)
		# Add in any needed triangles.
		for e1, e2 in combinations(self.edges, r=2):
			if e1 != e0 and e2 != e0:
				if e1.free_sides() > 0 and e2.free_sides() > 0:
					if len(set([e[0] for e in [e0, e1, e2]] + [e[1] for e in [e0, e1, e2]])) == 3:
						self.create_triangle(e0, e1, e2)
		self.unsaved_work = True
		self.redraw()
		self.build_equipped_triangulation()
		return self.edges[-1]
	
	def destroy_edge(self, edge=None):
		if edge is None: edge = self.edges[-1]
		if self.selected_object == edge: self.select_object(None)
		for drawn in edge.drawn:
			self.canvas.delete(drawn)
		for triangle in edge.in_triangles:
			self.destroy_triangle(triangle)
		self.destroy_edge_identification(edge)
		self.edges.remove(edge)
		self.unsaved_work = True
		self.build_equipped_triangulation()
		self.redraw()
	
	def create_triangle(self, e1, e2, e3):
		# Check that there are 3 edges.
		if len(set([e1, e2, e3])) != 3:
			return None
		
		# Check that this triangle doesn't already exist.
		if any([set(triangle.edges) == set([e1, e2, e3]) for triangle in self.triangles]):
			return None
		
		# Check that there are 3 vertices.
		corner_vertices = list(set(v for e in [e1, e2, e3] for v in e))
		if len(corner_vertices) != 3:
			return None
		
		# Check that there aren't any vertices inside the triangle.
		v0 = corner_vertices[2] - corner_vertices[0]
		v1 = corner_vertices[1] - corner_vertices[0]
		for vertex in self.vertices:
			if vertex not in corner_vertices:
				v2 = vertex - corner_vertices[0]
				
				dot00 = dot(v0, v0)
				dot01 = dot(v0, v1)
				dot02 = dot(v0, v2)
				dot11 = dot(v1, v1)
				dot12 = dot(v1, v2)
				
				invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01)
				u = (dot11 * dot02 - dot01 * dot12) * invDenom
				v = (dot00 * dot12 - dot01 * dot02) * invDenom
				
				if (u >= 0) and (v >= 0) and (u + v <= 1):
					return None
		
		self.triangles.append(curver.application.CanvasTriangle(self.canvas, [e1, e2, e3], self.options))
		
		self.unsaved_work = True
		self.redraw()
		self.build_equipped_triangulation()
		return self.triangles[-1]
	
	def destroy_triangle(self, triangle=None):
		if triangle is None: triangle = self.triangles[-1]
		self.canvas.delete(triangle.drawn)
		for edge in self.edges:
			if triangle in edge.in_triangles:
				edge.in_triangles.remove(triangle)
				self.destroy_edge_identification(edge)
		self.triangles.remove(triangle)
		self.unsaved_work = True
		self.redraw()
		self.build_equipped_triangulation()
	
	def create_edge_identification(self, e1, e2):
		if e1.equivalent_edge is not None or e2.equivalent_edge is not None:
			return None
		if e1.free_sides() != 1 or e2.free_sides() != 1:
			return None
		
		e1.equivalent_edge, e2.equivalent_edge = e2, e1
		# Now orient the edges so they match.
		[t1], [t2] = e1.in_triangles, e2.in_triangles
		[s1], [s2] = [i for i in range(3) if t1.edges[i] == e1], [i for i in range(3) if t2.edges[i] == e2]
		# Determine if the orientation of e1 (respectively e2) agrees with t1 (resp. t2).
		e1_agrees = e1[0] == t1[s1 + 1]
		e2_agrees = e2[0] == t2[s2 + 1]
		# We need one to agree and one to disagree - so if not then flip the orientation of e2.
		if e1_agrees == e2_agrees:
			e2.flip_orientation()
		
		# Change colour.
		new_colour = self.colour_picker.get_colour()
		e1.set_colour(new_colour)
		e2.set_colour(new_colour)
		self.unsaved_work = True
		self.build_equipped_triangulation()
	
	def destroy_edge_identification(self, edge):
		if edge.equivalent_edge is not None:
			other_edge = edge.equivalent_edge
			other_edge.set_colour()
			edge.set_colour()
			
			
			other_edge.set_colour(other_edge.default_colour)
			edge.set_colour(edge.default_colour)
			
			edge.equivalent_edge.equivalent_edge = None
			edge.equivalent_edge = None
			self.unsaved_work = True
		self.build_equipped_triangulation()
	
	def create_curve_component(self, vertices, multiplicity=1, smooth=False):
		self.curve_components.append(curver.application.CurveComponent(self.canvas, vertices, self.options, multiplicity, smooth))
		return self.curve_components[-1]
	
	def destory_curve_component(self, curve_component):
		if self.selected_object == curve_component: self.select_object(None)
		self.canvas.delete(curve_component.drawn)
		self.curve_components.remove(curve_component)
	
	def create_train_track_block(self, vertices, multiplicity=1, smooth=False):
		self.train_track_blocks.append(curver.application.TrainTrackBlock(self.canvas, vertices, self.options, multiplicity, smooth))
		return self.train_track_blocks[-1]
	
	def destroy_train_track_block(self, curve_component):
		self.canvas.delete(curve_component.drawn)
		self.train_track_blocks.remove(curve_component)
	
	def destroy_lamination(self):
		while self.curve_components != []:
			self.destory_curve_component(self.curve_components[-1])
		
		while self.train_track_blocks != []:
			self.destroy_train_track_block(self.train_track_blocks[-1])
		
		if self.is_complete():
			self.current_lamination = self.equipped_triangulation.triangulation.empty_lamination()
		
		self.select_object(None)
		self.redraw()
	
	
	######################################################################
	
	
	def set_edge_indices(self):
		# Assigns each edge an index in range(self.zeta).
		
		self.clear_edge_indices()
		self.zeta = 0
		for edge in self.edges:
			if edge.index == -1:
				self.zeta += 1
				edge.index = self.zeta-1
				if edge.equivalent_edge is not None:
					edge.equivalent_edge.index = edge.index
	
	def clear_edge_indices(self):
		self.zeta = 0
		for edge in self.edges:
			edge.index = -1
	
	def create_edge_labels(self):
		self.destroy_edge_labels()  # Remove existing labels.
		
		def accuracy_required(x):
			if isinstance(x, curver.IntegerType):
				return curver.kernel.height_int(x) + 1
			else:
				return x.height + x.log_degree
		def get_accurate(x, acc):
			if isinstance(x, curver.IntegerType):
				return curver.kernel.AlgebraicApproximation.from_int(x, acc)
			else:
				return x.algebraic_approximation(acc)
			
		
		# How to label the edge with given index.
		if self.options.label_edges == 'Index':
			labels = dict((index, index) for index in range(self.zeta))
		elif self.options.label_edges == 'Geometric':
			labels = dict((index, self.current_lamination(index)) for index in range(self.zeta))
		elif self.options.label_edges == 'Algebraic':
			labels = dict((index, self.current_lamination[index]) for index in range(self.zeta))
		elif self.options.label_edges == 'None':
			labels = dict((index, '') for index in range(self.zeta))
		else:
			raise ValueError()
		
		if self.options.projectivise and self.options.label_edges in ['Geometric', 'Algebraic']:
			accuracy_required = 2 * len(labels) * (max(accuracy_required(value) for value in labels.values()) + 1)
			labels = dict((index, get_accurate(labels[index], accuracy_required)) for index in labels)
			total = sum(labels[index] for index in range(self.zeta))
			if total != 0:  # There should probably be an else to this statement.
				# Note the "+ 0" to ensure that -0.0 appears as 0.0.
				labels = dict((index, round(float(labels[index] / total), 12) + 0) for index in range(self.zeta))
			else:
				labels = dict((index, round(float(labels[index]), 12) + 0) for index in range(self.zeta))
		
		for edge in self.edges:
			# We start by creating a nice background for the label. This ensures
			# that it is always readable, even when on top of a lamination.
			# To do this we first draw this label in a different colour with
			# slightly different offsets. This creates a nice 'bubble' effect
			# rather than having to draw a large bounding box.
			for offset in OFFSETS:
				self.canvas.create_text(
					[a+x for a, x in zip(edge.centre(), offset)],
					text=labels[edge.index],
					tag='edge_label',
					font=self.options.canvas_font,
					fill=DEFAULT_EDGE_LABEL_BG_COLOUR)
			
			self.canvas.create_text(edge.centre(),
				text=labels[edge.index],
				tag='edge_label',
				font=self.options.canvas_font,
				fill=DEFAULT_EDGE_LABEL_COLOUR)
	
	def destroy_edge_labels(self):
		self.canvas.delete('edge_label')
	
	def build_edge_labels(self):
		if self.is_complete():
			self.create_edge_labels()
		else:
			self.destroy_edge_labels()
	
	
	######################################################################
	
	
	def lamination_to_canvas(self, lamination):
		self.destroy_lamination()
		
		# We'll do everything with floats now because these are accurate enough for drawing to the screen with.
		vb = self.options.vertex_buffer  # We are going to use this a lot.
		a_weights = [float(x) for x in lamination]
		master = float(max(a_weights))
		if master == 0: master = float(1)
		
		for triangle in self.triangles:
			a_tri_weights = [a_weights[edge.index] for edge in triangle.edges]
			a_dual_weights = [(a_tri_weights[(j+1)%3] + a_tri_weights[(j+2)%3] - a_tri_weights[(j+0)%3]) / 2 for j in range(3)]
			for i in range(3):
				a = triangle[i-1] - triangle[i]
				b = triangle[i-2] - triangle[i]
				
				if lamination.weight() > MAX_DRAWABLE:
					if a_dual_weights[i] > 0.00001:  # Should be 0 but we have a floating point approximation.
						# We first do the edge to the left of the vertex.
						# Correction factor to take into account the weight on this edge.
						s_a = a_weights[triangle.edges[i-2].index] / master
						# The fractions of the distance of the two points on this edge.
						scale_a = vb * s_a + (1 - s_a) / 2
						scale_a2 = scale_a + (1 - 2*vb) * s_a * a_dual_weights[i] / (a_dual_weights[i] + a_dual_weights[i-1])
						
						# Now repeat for the other edge of the triangle.
						s_b = a_weights[triangle.edges[i-1].index] / master
						scale_b = vb * s_b + (1 - s_b) / 2
						scale_b2 = scale_b + (1 - 2*vb) * s_b * a_dual_weights[i] / (a_dual_weights[i] + a_dual_weights[i-2])
						
						S1, P1, Q1, E1 = curver.application.interpolate(triangle[i-1], triangle[i], triangle[i-2], scale_a, scale_b)
						S2, P2, Q2, E2 = curver.application.interpolate(triangle[i-1], triangle[i], triangle[i-2], scale_a2, scale_b2)
						self.create_train_track_block([S1, S1, P1, Q1, E1, E1, E2, E2, Q2, P2, S2, S2, S1, S1], smooth=True)
				else:
					# Also it is VERY slow (O(n) not O(log(n))).
					# Here we need the exact dual weights so we had better work them out.
					weights = [lamination(edge.index) for edge in triangle.edges]
					wa, wb = weights[i-2], weights[i-1]
					dual_weights = [(weights[j-2] + weights[j-1] - weights[j]) // 2 for j in range(3)]
					for j in range(int(dual_weights[i])):
						scale_a = 0.5 if wa == 1 else vb + (1 - 2*vb) * ((wa - 1) * (master - wa) + 2 * wa * j) / (2 * (wa - 1) * master)
						scale_b = 0.5 if wb == 1 else vb + (1 - 2*vb) * ((wb - 1) * (master - wb) + 2 * wb * j) / (2 * (wb - 1) * master)
						
						S, P, Q, E = curver.application.interpolate(triangle[i-1], triangle[i], triangle[i-2], scale_a, scale_b)
						if self.options.straight_laminations:
							self.create_curve_component([S, E])
						else:
							self.create_curve_component([S, P, Q, E], smooth=True)
		
		self.current_lamination = lamination
		self.create_edge_labels()
	
	
	######################################################################
	
	
	def parent_key_press(self, event):
		key = event.keysym
		focus = self.parent.focus_get()
		if key == 'F1':
			self.show_help()
		elif key == 'equal' or key == 'plus':
			self.zoom_in()
		elif key == 'minus' or key == 'underscore':
			self.zoom_centre(0.95)
		elif key == '0':
			self.zoom_to_drawing()
		elif key == 'Up':
			if focus == self.canvas:
				self.translate(0, 5)
		elif key == 'Down':
			if focus == self.canvas:
				self.translate(0, -5)
		elif key == 'Left':
			if focus == self.canvas:
				self.translate(5, 0)
		elif key == 'Right':
			if focus == self.canvas:
				self.translate(-5, 0)

def start(load_from=None):
	root = TK.Tk()
	root.title('curver')
	curver_application = CurverApplication(root)
	root.minsize(300, 300)
	root.geometry('700x500')
	root.wait_visibility(root)
	if load_from is not None: curver_application.load(load_from=load_from)
	# Set the icon.
	# Make sure to get the right path if we are in a cx_Freeze compiled executable.
	# See: http://cx-freeze.readthedocs.org/en/latest/faq.html#using-data-files
	datadir = os.path.dirname(sys.executable if getattr(sys, 'frozen', False) else __file__)
	icon_path = os.path.join(datadir, 'icon', 'icon.gif')
	img = TK.PhotoImage(file=icon_path)
	try:
		root.tk.call('wm', 'iconphoto', root._w, img)
	except TK.TclError:
		# Give up if we can't set the icon for some reason.
		# This seems to be a problem if you start curver within SnapPy.
		pass
	root.mainloop()
	return curver_application.output

if __name__ == '__main__':
	start()

