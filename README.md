# Raytracing project from scratch in c++


## Installation

You will find the project sources in the main folder, along the 3d models (.obj) and textures (.mtl + .bpm/.tga) in the folders next to it.

You will need to allow parallelism in your IDE for the project to work way faster than without.

## Details

This project has many features developed within which will be explained briefly here. To obtain more details about a part, you can check the code to understand it better or check on the Internet the theory behind each concept.

### Classic Geometry

The basic geometric features are coded in classes : Point/Vectors (3 coordinates in space), rays (1 point + 1 director vector), spheres (and later virtual objects), ...

The classic operations are also included, like scalar/vectorial operations, normalization, (square) norm, rotation in space (around a vector), multiplications by a scalar, ...

You will also find special methods within certain classes, like reflexion and refraction for the rays (linked to a normal vector).

### Objects intersections

One of the main operation needed for raytracing is to calcul the intersection between a ray and an object (a parallelepiped, a sphere, a cube, a triangle, ...).

Each object will then have a unique intersection method : second degree polynomial resolution for spheres, simple determinant for plans, barycentric coordinates for triangles, ...

### Light

The light on each pixel with be linked to the distance between the object intersected by the main ray and the light, but also by other factors like the scalar product between light/object direction and the ray direction.

### Shadows

A consequence of these intersections between objects and the light is the adding of basic shadows.
A way to have better defined shadows is to add diffusion to objects, via the Monte Carlo method to adapt the physics to computer calculs (which add a lot of time to render an image, consequence of the number of rays launched to a pixel of the image).
Finally, we can obtain soft shadows via the change of a single point light to a spherical light (working also with the Monte Carlo algorithm).

### Colors, mirror effect and transparency effect

Another feature enabled here is the coloring of objects, based on an unique albedo for each one of them (if wanted).

You will also be able to create mirror objects (from the reflexion basic physics behind) or even transparent ones (same as before, with basic refraction).

You can if you wish add these features to any objects (even mesh).

### Mesh loading and visualization

The next big thing about raytracing here is the loading of an already existing 3D model via a mesh (set of vertices and faces defining the model boundary). The main facing function to load a mesh (.off file) into our program is adapted from a code already existing, since not important to develop/understand for our project.

What is important is then to adapt the intersection method to a mesh (which is a set of triangles) : Intersecting a mesh is indeed intersecting a triangle within. It will takes naturally more time to calcul : A boundary box is set in place to reduce this time logarithmically (parallelepiped around other parallelepiped around ... around triangles, since parallelepipeds sintersections are faster to calcul than triangles').

### Textures loading

The final step in this project is the adding of textures for the models visualization. Same as the mesh loading function, the texture loading function is adapted from an already existing code.
The goal then is to adapt the color (or albedo) of a textured object to its texture, via the uv coordinates between the texture and the mesh.


## Exemples



