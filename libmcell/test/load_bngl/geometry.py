import mcell as m

# ---- Cube ----
Cube_vertex_list = [
    [-0.5, -0.5, -0.5], 
    [-0.5, -0.5, 0.5], 
    [-0.5, 0.5, -0.5], 
    [-0.5, 0.5, 0.5], 
    [0.5, -0.5, -0.5], 
    [0.5, -0.5, 0.5], 
    [0.5, 0.5, -0.5], 
    [0.5, 0.5, 0.5]
] # Cube_vertex_list

Cube_element_connections = [
    [1, 2, 0], 
    [3, 6, 2], 
    [7, 4, 6], 
    [5, 0, 4], 
    [6, 0, 2], 
    [3, 5, 7], 
    [1, 3, 2], 
    [3, 7, 6], 
    [7, 5, 4], 
    [5, 1, 0], 
    [6, 4, 0], 
    [3, 1, 5]
] # Cube_element_connections

Cube = m.GeometryObject(
    name = 'Cube',
    vertex_list = Cube_vertex_list,
    element_connections = Cube_element_connections,
    surface_regions = []
)
# ^^^^ Cube ^^^^


