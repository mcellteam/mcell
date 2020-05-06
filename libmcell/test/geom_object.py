#!/usr/bin/env python3

import mcell as m


tetrahedron_vert_list = [
    [  0.00,  0.00,  0.02 ],
    [  0.02,  0.00, -0.01 ],
    [ -0.01,  0.02, -0.01 ],
    [ -0.01, -0.02, -0.01 ]
]

tetrahedron_face_list = [
    [ 0, 1, 2 ],
    [ 0, 2, 3 ],
    [ 0, 3, 1 ],
    [ 1, 3, 2 ]
]    

tetrahedron_name_surface_list = [
    0, 1
] 

sr = m.SurfaceRegion('name', tetrahedron_name_surface_list)

print(sr)

obj = m.GeometryObject(
    'tetrahedron', 
    tetrahedron_vert_list, tetrahedron_face_list, 
    [sr]
)

print(obj)

instantiation = m.InstantiationData()

instantiation.instantiate_geometry_object(obj)


