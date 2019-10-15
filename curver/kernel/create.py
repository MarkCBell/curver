
''' A module for providing creation functions for all of the moves within curver. '''

import curver

def link(move, inverse_move):
    ''' Link two moves to be each others inverses and return the first. '''
    move._inverse = inverse_move
    inverse_move._inverse = move
    return move

def isometry(source_triangulation, target_triangulation, label_map):
    ''' Create an isometry. '''
    inverse_label_map = dict((label_map[label], label) for label in label_map)
    return link(
        curver.kernel.Isometry(source_triangulation, target_triangulation, label_map),
        curver.kernel.Isometry(target_triangulation, source_triangulation, inverse_label_map)
        )

def edgeflip(source_triangulation, target_triangulation, edge):
    ''' Create an edgeflip. '''
    return link(
        curver.kernel.EdgeFlip(source_triangulation, target_triangulation, edge),
        curver.kernel.EdgeFlip(target_triangulation, source_triangulation, ~edge)
        )

def multiedgeflip(source_triangulation, target_triangulation, edges):
    ''' Create a multiedgeflip. '''
    return link(
        curver.kernel.MultiEdgeFlip(source_triangulation, target_triangulation, edges),
        curver.kernel.MultiEdgeFlip(target_triangulation, source_triangulation, [~edge for edge in edges])
        )

def twist(curve, power):
    ''' Create a twist. '''
    return link(
        curver.kernel.Twist(curve, power),
        curver.kernel.Twist(curve, -power)
        )

def halftwist(curve, power):
    ''' Create a halftwist. '''
    return link(
        curver.kernel.HalfTwist(curve, power),
        curver.kernel.HalfTwist(curve, -power)
        )

def crush(source_triangulation, target_triangulation, curve, matrix):
    ''' Create a crush. '''
    return link(
        curver.kernel.Crush(source_triangulation, target_triangulation, curve),
        curver.kernel.Lift(target_triangulation, source_triangulation, matrix)  # pylint: disable=arguments-out-of-order
        )

def lineartransformation(source_triangulation, target_triangulation, matrix, inverse_matrix):
    ''' Create a linear transformation. '''
    return link(
        curver.kernel.LinearTransformation(source_triangulation, target_triangulation, matrix),
        curver.kernel.LinearTransformation(target_triangulation, source_triangulation, inverse_matrix)
        )

