# GJKEPA
Pure C# implementation of GJK and EPA. Fast and stable, with a minimum of external links. 

Provides the determination of the fact of intersection between two polyhedra.
Upon intersection, it returns the average contact point, the collision exit vector, and the penetration depth.
The implementation of GJK was taken and adapted from the article https://winter.dev/articles/gjk-algorithm.
The implementation of EPA was taken and adapted from the article https://winter.dev/articles/epa-algorithm.
To determine the contact point, the barycentric coordinates of the projection of the origin were used.
