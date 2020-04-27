#!/usr/bin/env python3

import mcell as m


tetrahedron_vert_list = [
    [  0.00,  0.00,  0.02, 0.5 ],
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

obj = m.GeometryObject(
    'tetrahedron', 
    tetrahedron_vert_list, tetrahedron_face_list
)

print(obj)
