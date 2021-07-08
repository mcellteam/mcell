"""
This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to [http://unlicense.org] 
"""

# Bring in Blender's Python API
import bpy



def add_mesh(verts, faces, name):
    mesh = bpy.data.meshes.new(name)
    mesh.from_pydata(verts, [], faces)
    mesh.update()

    obj = bpy.data.objects.new("Obj_" + name, mesh)
    scene = bpy.context.scene
    scene.objects.link(obj)
    return obj


def add_point(pos, name):
    point_to_triangle = [
      # insert position of the molecule here
      pos,
      ( pos[0], pos[1] + 0.01, pos[2] + 10),
      ( pos[0], pos[1] + 0.01, pos[2] - 0.01),
    ] 

    point_faces = [  
      (0, 1, 2)
    ]
        
    pt = bpy.data.meshes.new("molecule")
    pt.from_pydata(point_to_triangle, [], point_faces)
    pt.update()

    obj = bpy.data.objects.new("Obj_" + name, pt)
    scene = bpy.context.scene
    scene.objects.link(obj)
    return obj
    
# add generated code below
mesh_verts = [
  (0, 0, 2), #0
  (-1, -2, -1), #1
  (2, 0, -1), #2
  (0, 0, 1), #3
  (-1, -2, -1), #4
  (2, 0, -1), #5
]
mesh_faces = [
  (0, 1, 2, ),
  (3, 4, 5, ),
  (0, 1, 3, ),
  (3, 1, 4, ),
  (1, 2, 4, ),
  (4, 2, 5, ),
  (2, 0, 5, ),
  (5, 0, 3, ),
]
add_mesh(mesh_verts, mesh_faces, "my_mesh")
mol0 = (0.583606, 0.255219, 0.63674)
add_point(mol0, "molecule0")
mol1 = (0.413867, 0.175317, -0.334765)
add_point(mol1, "molecule1")
mol2 = (0.195182, -0.446152, 0.66672)
add_point(mol2, "molecule2")
mol3 = (-0.47739, 0.532683, -0.000138335)
add_point(mol3, "molecule3")

