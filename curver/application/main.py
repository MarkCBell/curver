
''' The main window of the curver GUI application. '''

import curver
import curver.application

import os
import sys
from math import sin, cos, pi, ceil, sqrt
from random import random
from colorsys import hls_to_rgb
from collections import namedtuple

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

Showable = namedtuple('Showable', ['name', 'item'])

def dot(a, b):
    return a[0] * b[0] + a[1] * b[1]

PHI = (1 + sqrt(5)) / 2

def get_colours(num_colours):
    def colour(state):
        hue = (state / PHI) % 1.
        lightness = (50 + random() * 10)/100.
        saturation = (90 + random() * 10)/100.
        r, g, b = hls_to_rgb(hue, lightness, saturation)
        return '#%02x%02x%02x' % (int(r * 255), int(g * 255), int(b * 255))
    
    return [colour(i) for i in range(num_colours)]

class Drawing(object):
    def __init__(self, note, canvas, item, options):
        self.note = note
        self.canvas = canvas
        self.options = options
        self.vertices = []
        self.edges = []
        self.triangles = []
        self.curve_components = []
        self.item = item
        self.lamination = None
        
        if isinstance(item, curver.kernel.MappingClassGroup):
            self.draw_lamination(item.triangulation.empty_lamination())
        if isinstance(item, curver.kernel.Triangulation):
            self.draw_lamination(item.empty_lamination())
        elif isinstance(item, curver.kernel.Lamination):
            self.draw_lamination(item)
        elif isinstance(item, curver.kernel.MappingClass):
            self.draw_lamination(item(item.source_triangulation.as_lamination()))
    
    def get_size(self):
        # Return the size of the canvas.
        # We would like to do:
        # return int(self.canvas.winfo_width()), int(self.canvas.winfo_height())
        # But the canvas may not have been rendered yet so instead we use the note:
        
        return (self.note.winfo_width() - 54, self.note.winfo_height() - 62)
    
    ######################################################################
    
    def create_vertex(self, point):
        assert isinstance(point, curver.application.Vector2)
        self.vertices.append(curver.application.CanvasVertex(self.canvas, point, self.options))
        return self.vertices[-1]
    
    def create_edge(self, v1, v2, label, colour, create_inverse=False):
        if create_inverse: self.create_edge(v2, v1, ~label, colour)
        self.edges.append(curver.application.CanvasEdge(self.canvas, [v1, v2], label, colour, self.options))
        return self.edges[-1]
    
    def create_triangle(self, e1, e2, e3):
        self.triangles.append(curver.application.CanvasTriangle(self.canvas, [e1, e2, e3], self.options))
        return self.triangles[-1]
    
    def create_curve_component(self, vertices, thin=True, smooth=False):
        assert all(isinstance(vertex, curver.application.Vector2) for vertex in vertices)
        self.curve_components.append(curver.application.CurveComponent(self.canvas, vertices, self.options, thin, smooth))
        return self.curve_components[-1]
    
    ######################################################################
    
    def translate(self, dx, dy):
        dv = curver.application.Vector2(dx, dy)
        for vertex in self.vertices:
            vertex.vector = vertex.vector + dv
        
        for curve_component in self.curve_components:
            for i in range(len(curve_component.vertices)):
                curve_component.vertices[i] = curve_component.vertices[i] + dv
        
        self.canvas.move('all', dx, dy)
    
    def zoom(self, scale):
        for vertex in self.vertices:
            vertex.vector = scale*vertex.vector
            vertex.update()
        for edge in self.edges:
            edge.update()
        for triangle in self.triangles:
            triangle.update()
        for curve_component in self.curve_components:
            for i in range(len(curve_component.vertices)):
                curve_component.vertices[i] = scale * curve_component.vertices[i]
            curve_component.update()
        self.redraw()
    
    def zoom_centre(self, scale):
        cw, ch = self.get_size()
        self.translate(-cw / 2, -ch / 2)
        self.zoom(scale)
        self.translate(cw / 2, ch / 2)
    
    def zoom_in(self):
        self.zoom_centre(1.05)
    
    def zoom_out(self):
        self.zoom_centre(0.95)
    
    def zoom_to_drawing(self):
        box = self.canvas.bbox('all')
        if box is not None:
            x0, y0, x1, y1 = box
            cw, ch = self.get_size()
            cr = min(cw, ch)
            
            w, h = x1 - x0, y1 - y0
            r = max(w, h)
            
            self.translate(-x0 - w / 2, -y0 - h / 2)
            self.zoom(self.options.zoom_fraction * float(cr) / r)
            self.translate(cw / 2, ch / 2)
    
    ######################################################################
    
    def draw_triangulation(self, triangulation):
        # Get a dual tree.
        dual_tree = triangulation.dual_tree()
        colours = dict((index, None) for index in triangulation.indices)
        outside = [index for index in triangulation.indices if not dual_tree[index]]
        for index, colour in zip(outside, get_colours(len(outside))):
            colours[index] = colour
        components = triangulation.components()
        num_components = len(components)
        cw, ch = self.get_size()
        
        # We will layout the components in a p x q grid.
        # Aim to maximise r := min(cw / p, ch / q) subject to pq >= num_components.
        # There is probably a closed formula for the optimal value of p (and so q).
        p = max(range(1, num_components+1), key=lambda p: min(cw / p, ch / ceil(float(num_components) / p)))
        q = int(ceil(float(num_components) / p))
        
        r = min(cw / p, ch / q) * (1 + self.options.zoom_fraction) / 4
        dx = cw / p
        dy = ch / q
        
        num_used_vertices = 0
        for index, component in enumerate(components):
            n = len(component) // 3  # Number of triangles.
            ngon = n + 2
            
            # Create the vertices.
            for i in range(ngon):
                self.create_vertex(
                    curver.application.Vector2(
                        dx * (index % p) + dx / 2 + r * sin(2 * pi * (i + 0.5) / ngon),
                        dy * int(index / p) + dy / 2 + r * cos(2 * pi * (i + 0.5) / ngon)
                    ))
            
            def num_descendants(edge_label):
                ''' Return the number of triangles that can be reached in the dual tree starting at the given edge_label. '''
                
                corner = triangulation.corner_lookup[edge_label]
                left = (1 + sum(num_descendants(~(corner.labels[2])))) if dual_tree[corner.indices[2]] else 0
                right = (1 + sum(num_descendants(~(corner.labels[1])))) if dual_tree[corner.indices[1]] else 0
                
                return left, right
            
            initial_edge_index = min(edge.index for edge in component if not dual_tree[edge.index])
            to_extend = [(num_used_vertices, num_used_vertices+1, initial_edge_index)]
            # Hmmm, need to be more careful here to ensure that we correctly orient the edges.
            self.create_edge(self.vertices[num_used_vertices+1], self.vertices[num_used_vertices+0], initial_edge_index, colours[initial_edge_index])
            while to_extend:
                source_vertex, target_vertex, label = to_extend.pop()
                left, right = num_descendants(label)
                far_vertex = target_vertex + left + 1
                corner = triangulation.corner_lookup[label]
                
                if corner[2].sign() == +1:
                    self.create_edge(self.vertices[far_vertex], self.vertices[target_vertex], corner[2].label, colours[corner[2].index], left > 0)
                else:
                    self.create_edge(self.vertices[target_vertex], self.vertices[far_vertex], corner[2].label, colours[corner[2].index], left > 0)
                if corner[1].sign() == +1:
                    self.create_edge(self.vertices[source_vertex], self.vertices[far_vertex], corner[1].label, colours[corner[1].index], right > 0)
                else:
                    self.create_edge(self.vertices[far_vertex], self.vertices[source_vertex], corner[1].label, colours[corner[1].index], right > 0)
                
                if left > 0:
                    to_extend.append((far_vertex, target_vertex, ~(corner[2].label)))
                
                if right > 0:
                    to_extend.append((source_vertex, far_vertex, ~(corner[1].label)))
            num_used_vertices = len(self.vertices)
        
        self.edges = sorted(self.edges, key=lambda e: ((0 if e.label >= 0 else 1), e.label))  # So now self.edges[label] is the edge with label label.
        
        for triangle in triangulation:
            self.create_triangle(self.edges[triangle[0].label], self.edges[triangle[1].label], self.edges[triangle[2].label])
    
    def draw_lamination(self, lamination):
        self.draw_triangulation(lamination.triangulation)  # This starts with self.initialise().
        
        vb = self.options.vertex_buffer  # We are going to use this a lot.
        master = float(max(abs(weight) for weight in lamination))
        
        if master > MAX_DRAWABLE:
            for triangle in self.triangles:
                weights = [max(lamination(edge.label), 0) for edge in triangle.edges]
                dual_weights = [lamination.dual_weight(edge.label) for edge in triangle.edges]
                parallel_arcs = [max(-lamination(edge.label), 0) for edge in triangle.edges]
                parallel_weights = [weight // 2 + (weight % 2 if edge.label >= 0 else 0) for edge, weight in zip(triangle.edges, parallel_arcs)]
                for i in range(3):
                    # Dual arcs.
                    if dual_weights[i] > 0:
                        # We first do the edge to the left of the vertex.
                        # Correction factor to take into account the weight on this edge.
                        s_a = (1 - 2*vb) * weights[i-2] / master
                        # The fractions of the distance of the two points on this edge.
                        scale_a = (1 - s_a) / 2
                        scale_a2 = scale_a + s_a * dual_weights[i] / weights[i-2]
                        
                        # Now repeat for the other edge of the triangle.
                        s_b = (1 - 2*vb) * weights[i-1] / master
                        scale_b = (1 - s_b) / 2
                        scale_b2 = scale_b + s_b * dual_weights[i] / weights[i-1]
                        
                        S1, P1, Q1, E1 = curver.application.interpolate(triangle[i-1].vector, triangle[i].vector, triangle[i-2].vector, scale_a, scale_b)
                        S2, P2, Q2, E2 = curver.application.interpolate(triangle[i-1].vector, triangle[i].vector, triangle[i-2].vector, scale_a2, scale_b2)
                        self.create_curve_component([S1, S1, P1, Q1, E1, E1, E2, E2, Q2, P2, S2, S2, S1, S1], thin=False)
                    elif dual_weights[i] < 0:  # Terminal arc.
                        s_0 = (1 - 2*vb) * weights[i] / master
                        
                        scale_a = (1 - s_0) / 2 + s_0 * dual_weights[i-2] / weights[i]
                        scale_a2 = scale_a + s_0 * (-dual_weights[i]) / weights[i]
                        
                        S1, P1, Q1, E1 = curver.application.interpolate(triangle[i-2].vector, triangle[i-1].vector, triangle[i].vector, scale_a, 1.0)
                        S2, P2, Q2, E2 = curver.application.interpolate(triangle[i-2].vector, triangle[i-1].vector, triangle[i].vector, scale_a2, 1.0)
                        self.create_curve_component([S1, S1, P1, E1, E1, P2, S2, S2, S1, S1], thin=False)
                    else:  # dual_weights[i] == 0:  # Nothing to draw.
                        pass
                    
                    # Parallel arcs.
                    if parallel_weights[i]:
                        S, O, E = triangle[i-2].vector, triangle[i].vector, triangle[i-1].vector
                        M = (S + E) / 2.0
                        centroid = (S + O + E) / 3.0
                        P = M + 1.0 * (centroid - M) * (parallel_weights[i]+1) / master
                        self.create_curve_component([S, S, P, E, E, S, S], thin=False)
        else:  # Draw everything. Caution, this is is VERY slow (O(n) not O(log(n))) so we only do it when the weight is low.
            for triangle in self.triangles:
                weights = [max(lamination(edge.label), 0) for edge in triangle.edges]
                dual_weights = [lamination.dual_weight(edge.label) for edge in triangle.edges]
                parallel_arcs = [max(-lamination(edge.label), 0) for edge in triangle.edges]
                parallel_weights = [weight // 2 + (weight % 2 if edge.label >= 0 else 0) for edge, weight in zip(triangle.edges, parallel_arcs)]
                for i in range(3):  # Dual arcs:
                    if dual_weights[i] > 0:
                        s_a = (1 - 2*vb) * weights[i-2] / master
                        s_b = (1 - 2*vb) * weights[i-1] / master
                        for j in range(dual_weights[i]):
                            scale_a = 0.5 if weights[i-2] == 1 else (1 - s_a) / 2 + s_a * j / (weights[i-2] - 1)
                            scale_b = 0.5 if weights[i-1] == 1 else (1 - s_b) / 2 + s_b * j / (weights[i-1] - 1)
                            
                            S, P, Q, E = curver.application.interpolate(triangle[i-1].vector, triangle[i].vector, triangle[i-2].vector, scale_a, scale_b)
                            self.create_curve_component([S, P, Q, E])
                    elif dual_weights[i] < 0:  # Terminal arc.
                        s_0 = (1 - 2*vb) * weights[i] / master
                        for j in range(-dual_weights[i]):
                            scale_a = 0.5 if weights[i] == 1 else (1 - s_0) / 2 + s_0 * dual_weights[i-1] / (weights[i] - 1) + s_0 * j / (weights[i] - 1)
                            
                            S, P, Q, E = curver.application.interpolate(triangle[i-2].vector, triangle[i-1].vector, triangle[i].vector, scale_a, 1.0)
                            self.create_curve_component([S, P, E])
                    else:  # dual_weights[i] == 0:  # Nothing to draw.
                        pass
                    
                    # Parallel arcs:
                    S, O, E = triangle[i-2].vector, triangle[i].vector, triangle[i-1].vector
                    M = (S + E) / 2.0
                    centroid = (S + O + E) / 3.0
                    for j in range(parallel_weights[i]):
                        P = M + float(j+1) / (parallel_weights[i] + 1) * (centroid - M) * (parallel_weights[i] + 1) / master
                        self.create_curve_component([S, P, E])
        
        self.lamination = lamination
        self.zoom_to_drawing()  # Recheck.
        self.redraw()  # Install options.
    
    def redraw(self):
        self.create_edge_labels()
        
        self.canvas.itemconfig('line', width=self.options.line_size)
        self.canvas.itemconfig('curve', width=self.options.line_size)
        # Only put arrows on the start so arrow heads appear in the middle.
        for edge in self.edges:
            if self.options.show_orientations and (edge.label >= 0 or self.edges[~edge.label].vertices[::-1] != edge.vertices):
                self.canvas.itemconfig(edge.drawn[0], arrow='last', arrowshape=self.options.arrow_shape)
            else:
                self.canvas.itemconfig(edge.drawn[0], arrow='')
        
        for edge in self.edges:
            edge.hide(not self.options.show_internals and self.edges[~edge.label].vertices[::-1] == edge.vertices)
        self.canvas.tag_raise('polygon')
        self.canvas.tag_raise('line')
        self.canvas.tag_raise('oval')
        self.canvas.tag_raise('curve')
        self.canvas.tag_raise('label')
        self.canvas.itemconfig('curve', smooth=self.options.smooth)
    
    ######################################################################
    
    def destroy_edge_labels(self):
        self.canvas.delete('label')
    
    def create_edge_labels(self):
        self.destroy_edge_labels()  # Remove existing labels.
        
        # How to label each edge.
        if self.options.label_edges == curver.application.options.LABEL_EDGES_INDEX:
            labels = dict((edge.label, edge.index) for edge in self.lamination.triangulation.edges)
        elif self.options.label_edges == curver.application.options.LABEL_EDGES_GEOMETRIC:
            labels = dict((edge.label, self.lamination(edge)) for edge in self.lamination.triangulation.edges)
        elif self.options.label_edges == curver.application.options.LABEL_EDGES_GEOMETRIC_PROJ:
            labels = dict((edge.label, self.lamination(edge)) for edge in self.lamination.triangulation.edges)
            total = float(sum(labels.values()))
            if total != 0:
                # Note the "+ 0" to ensure that -0.0 appears as 0.0.
                labels = dict((label, round(labels[label] / total, 12) + 0) for label in labels)
            else:
                labels = dict((label, 0.0) for label in labels)
        elif self.options.label_edges == curver.application.options.LABEL_EDGES_NONE:
            labels = dict((edge.label, '') for edge in self.lamination.triangulation.edges)
        else:
            raise ValueError()
        
        for edge in self.edges:
            if edge.label >= 0 or self.edges[~edge.label].vertices[::-1] != edge.vertices:
                # We start by creating a nice background for the label. This ensures
                # that it is always readable, even when on top of a lamination.
                # To do this we first draw this label in a different colour with
                # slightly different offsets. This creates a nice 'bubble' effect
                # rather than having to draw a large bounding box.
                for offset in OFFSETS:
                    self.canvas.create_text(
                        [a+x for a, x in zip(edge.centre().to_tuple(), offset)],
                        text=labels[edge.label],
                        tag='label',
                        font=self.options.canvas_font,
                        fill=DEFAULT_EDGE_LABEL_BG_COLOUR
                        )
                
                self.canvas.create_text(
                    edge.centre().to_tuple(),
                    text=labels[edge.label],
                    tag='label',
                    font=self.options.canvas_font,
                    fill=DEFAULT_EDGE_LABEL_COLOUR
                    )


class CurverApplication(object):
    def __init__(self, parent, items):
        self.parent = parent
        self.options = curver.application.Options(self)
        
        s = TTK.Style()
        s.configure('TNotebook', tabposition='n')
        
        self.note = TTK.Notebook(self.parent)
        self.note.pack(fill='both', expand=True, padx=2, pady=2)
        self.note.enable_traversal()
        
        if isinstance(items, (list, tuple)) and len(items) == 1 and isinstance(items[0], (list, dict)):  # [list / dict].
            items = items[0]
        
        if isinstance(items, dict):  # dict.
            items = [Showable(key, item) for key, item in sorted(dict(items).items(), key=lambda x: x[0])]
        try:  # list of pairs.
            items = [Showable(key, item) for key, item in items]
        except (TypeError, ValueError):  # list.
            items = [Showable(name, item) for name, item in curver.kernel.utilities.name_objects(items)]
        self.items = items
        
        self.canvases = [TK.Canvas(self.note, bg='#dcecff') for i in range(len(self.items))]
        for canvas, item in zip(self.canvases, self.items):
            self.note.add(canvas, text=item.name)
        
        self.parent.wait_visibility(self.note)  # Make sure we get the right sizes of things.
        self.drawings = [Drawing(self.note, canvas, item.item, self.options) for canvas, item in zip(self.canvases, self.items)]
        
        # Create the menus.
        # Make sure to start the Lamination and Mapping class menus disabled.
        self.menubar = TK.Menu(self.parent)
        
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
        self.edgelabelmenu.add_radiobutton(label=curver.application.options.LABEL_EDGES_GEOMETRIC_PROJ, var=self.options.label_edges_var)
        
        self.zoommenu = TK.Menu(self.menubar, tearoff=0)
        self.zoommenu.add_command(label='Zoom in', command=self.zoom_in, accelerator='+')
        self.zoommenu.add_command(label='Zoom out', command=self.zoom_out, accelerator='-')
        self.zoommenu.add_command(label='Zoom to drawing', command=self.zoom_to_drawing, accelerator='0')
        
        self.settingsmenu.add_cascade(label='Sizes', menu=self.sizemenu)
        self.settingsmenu.add_cascade(label='Edge label', menu=self.edgelabelmenu)
        self.settingsmenu.add_cascade(label='Zoom', menu=self.zoommenu)
        self.settingsmenu.add_checkbutton(label='Smooth', var=self.options.smooth_var)
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
    
    def current_drawing(self):
        return self.drawings[self.note.index('current')]
    
    def redraw(self):
        for drawing in self.drawings:
            drawing.redraw()
    
    def zoom_in(self):
        self.current_drawing().zoom_in()
    
    def zoom_out(self):
        self.current_drawing().zoom_out()
    
    def zoom_to_drawing(self):
        self.current_drawing().zoom_to_drawing()
    
    def translate(self, dx, dy):
        self.current_drawing().translate(dx, dy)
    
    def parent_key_press(self, event):
        key = event.keysym
        if key == 'F1':
            self.show_help()
        elif key == 'equal' or key == 'plus':
            self.zoom_in()
        elif key == 'minus' or key == 'underscore':
            self.zoom_out()
        elif key == '0':
            self.zoom_to_drawing()
        elif key == 'Up':
            self.translate(0, 5)
        elif key == 'Down':
            self.translate(0, -5)
        elif key == 'Left':
            self.translate(5, 0)
        elif key == 'Right':
            self.translate(-5, 0)
    
    def export_image(self):
        path = tkFileDialog.asksaveasfilename(defaultextension='.ps', filetypes=[('postscript files', '.ps'), ('all files', '.*')], title='Export Image')
        if path:
            try:
                self.current_drawing().canvas.postscript(file=path, colormode='color')
            except IOError:
                tkMessageBox.showwarning('Export Error', 'Could not open: %s' % path)
    
    def quit(self):
        # Apparantly there are some problems with comboboxes, see:
        #  http://stackoverflow.com/questions/15448914/python-tkinter-ttk-combobox-throws-exception-on-quit
        self.parent.eval('::ttk::CancelRepeat')
        self.parent.destroy()
        self.parent.quit()
    
    def show_help(self):
        pass
        # curver.doc.open_documentation()
    
    def show_about(self):
        tkMessageBox.showinfo('About', 'curver (Version %s).\nCopyright (c) Mark Bell 2017.' % curver.__version__)

def start(*items):
    root = TK.Tk()
    root.title('curver')
    root.minsize(300, 300)
    root.geometry('700x500')
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
    CurverApplication(root, items)
    root.mainloop()

if __name__ == '__main__':
    start()

