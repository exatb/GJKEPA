using System;
using System.Collections.Generic;
using System.Numerics;
using System.Linq;

/// <summary>
/// Обеспечивает определение факта пересечения двух многогранников. 
/// При пересечении возвращает среднюю точку контакта, вектор выхода из коллизии и глубину проникновения. 
/// Реализация GJK взята и доработана из статьи https://winter.dev/articles/gjk-algorithm.
/// Реализация EPA взята и доработана из статьи https://winter.dev/articles/epa-algorithm.
/// Для опрделения точки контакта вычисляются барицентрические координаты проекции начала координат  
/// на возвращенный EPA ближайший к началу координат треугольник и переводятся в координаты исходных многогранников.
/// </summary>

public class GJK_EPA_BСP
{
    static public int gjk_max = 32; //Maximum gjk iterations
    static public int epa_max = 32; //Maximum epa iterations

    private const float Epsilon = 1e-3f; // Малое значение для численной стабильности EPA

#if GJKEPA_DEBUG
    public delegate void DrawPointDelegate(Vector3 v, UInt32 color);
    public delegate void DrawVectorDelegate(Vector3 begin, Vector3 end, UInt32 color);
    public delegate void DrawLineDelegate(Vector3 begin, Vector3 end, UInt32 color);
    public delegate void LogDelegate(string s);

    //Debugging functions will run without checking for null!
    static public DrawPointDelegate DrawPoint = null;
    static public DrawVectorDelegate DrawVector = null;
    static public DrawLineDelegate DrawLine = null;
    static public LogDelegate Log = null;
#endif  

    public static bool CheckIntersection(Vector3[] polyhedron1, Vector3[] polyhedron2, out Vector3 contactPoint, out float penetrationDepth, out Vector3 contactNormal)
    {
        contactPoint = new Vector3(0, 0, 0);
        penetrationDepth = 0;
        contactNormal = new Vector3(0, 0, 0);

        // Инициализация GJK
        Simplex simplex = new Simplex();
        Vector3 direction = Vector3.Subtract(polyhedron1[0], polyhedron2[0]);

        SupportPoint support = CalculateSupportPoint(polyhedron1, polyhedron2, direction);

        simplex.AddSupportPoint(support);

        direction = -support.Point;

        int gjk_cnt = 0;

        while (gjk_cnt < gjk_max)
        {
            support = CalculateSupportPoint(polyhedron1, polyhedron2, direction);

            if (Vector3.Dot(support.Point, direction) <= 0)
            {
                return false; // Нет пересечения
            }

            simplex.AddSupportPoint(support);

            if (simplex.ContainsOrigin(ref direction))
            {
                // Обнаружено пересечение, продолжаем с EPA
#if GJKEPA_DEBUG                
                //simplex.DebugDraw();
#endif

                int epa_cnt = EPA(polyhedron1, polyhedron2, simplex, out contactNormal, out penetrationDepth, out contactPoint);

#if GJKEPA_DEBUG                
                DrawVector(contactPoint, contactNormal * penetrationDepth, 4278222848u); //Green out vector
                DrawPoint(contactPoint, 4294901760u); //Red pont
                Log(" gjk=" + gjk_cnt.ToString() + "  , epa=" + epa_cnt.ToString());
#endif
                return true;
            }

            gjk_cnt++;
        }
#if GJKEPA_DEBUG
        Log(" gjk=" + gjk_cnt.ToString());
#endif
        return false;
    }

    public struct SupportPoint
    {
        public Vector3 Point; // Точка поддержки
        public int Index1;    // Индекс в первом многограннике
        public int Index2;    // Индекс во втором многограннике

        public SupportPoint(Vector3 point, int index1, int index2)
        {
            Point = point;
            Index1 = index1;
            Index2 = index2;
        }
    }

    private static SupportPoint CalculateSupportPoint(Vector3[] polyhedron1, Vector3[] polyhedron2, Vector3 direction)
    {
        int index1, index2;
        Vector3 support1 = GetSupportPoint(polyhedron1, direction, out index1);
        Vector3 support2 = GetSupportPoint(polyhedron2, -direction, out index2);

        return new SupportPoint(Vector3.Subtract(support1, support2), index1, index2);
    }

    private static Vector3 GetSupportPoint(Vector3[] polyhedron, Vector3 direction, out int index)
    {
        double maxDot = double.MinValue;
        Vector3 support = new Vector3(0, 0, 0);
        index = -1;

        for (int i = 0; i < polyhedron.Length; i++)
        {
            double dot = Vector3.Dot(polyhedron[i], direction);
            if (dot > maxDot)
            {
                maxDot = dot;
                support = polyhedron[i];
                index = i;
            }
        }

        return support;
    }

    private class Simplex
    {
        private List<SupportPoint> points;

        public List<SupportPoint> Vertices
        {
            get { return points; }
        }

        public Simplex()
        {
            points = new List<SupportPoint>();
        }

#if GJKEPA_DEBUG                
        public void DebugDraw()
        {
            UInt32 clr = 4278222848u; //Green color
            float mul = 1.05f;
            switch (points.Count)
            {
                case 1:
                    DrawPoint(points[0].Point, clr);
                    break;
                case 2:
                    DrawLine(points[0].Point * mul, points[1].Point * mul, clr);
                    break;
                case 3:
                    DrawLine(points[0].Point * mul, points[1].Point * mul, clr);
                    DrawLine(points[1].Point * mul, points[2].Point * mul, clr);
                    DrawLine(points[2].Point * mul, points[0].Point * mul, clr);
                    break;
                case 4:
                    DrawLine(points[0].Point * mul, points[1].Point * mul, clr);
                    DrawLine(points[0].Point * mul, points[2].Point * mul, clr);
                    DrawLine(points[0].Point * mul, points[3].Point * mul, clr);
                    DrawLine(points[1].Point * mul, points[3].Point * mul, clr);
                    DrawLine(points[1].Point * mul, points[3].Point * mul, clr);
                    DrawLine(points[2].Point * mul, points[3].Point * mul, clr);
                    DrawLine(points[1].Point * mul, points[2].Point * mul, clr);
                    break;
                default:
                    return;
            }
        }
#endif

        public void AddSupportPoint(SupportPoint supportPoint)
        {
            points.Insert(0, supportPoint);
        }

        public bool ContainsOrigin(ref Vector3 direction)
        {
            if (points.Count == 4)
            {
                return ContainsOriginTetrahedron(ref direction);
            }
            else if (points.Count == 3)
            {
                return ContainsOriginTriangle(ref direction);
            }
            else if (points.Count == 2)
            {
                return ContainsOriginLine(ref direction);
            }

            return false;
        }

        private bool ContainsOriginTetrahedron(ref Vector3 direction)
        {
            SupportPoint a = points[0];
            SupportPoint b = points[1];
            SupportPoint c = points[2];
            SupportPoint d = points[3];

            Vector3 ab = b.Point - a.Point;
            Vector3 ac = c.Point - a.Point;
            Vector3 ad = d.Point - a.Point;
            Vector3 ao = -a.Point;

            Vector3 abc = Vector3.Cross(ab, ac);
            Vector3 acd = Vector3.Cross(ac, ad);
            Vector3 adb = Vector3.Cross(ad, ab);

            if (Vector3.Dot(abc, ao) > 0)
            {
                points.Clear();
                points.Add(a);
                points.Add(b);
                points.Add(c);
                //return Triangle(points = { a, b, c }, direction);
                return ContainsOriginTriangle(ref direction);
            }

            if (Vector3.Dot(acd, ao) > 0)
            {
                points.Clear();
                points.Add(a);
                points.Add(c);
                points.Add(d);
                //return Triangle(points = { a, c, d }, direction);
                return ContainsOriginTriangle(ref direction);
            }

            if (Vector3.Dot(adb, ao) > 0)
            {
                points.Clear();
                points.Add(a);
                points.Add(d);
                points.Add(b);
                //return Triangle(points = { a, d, b }, direction);
                return ContainsOriginTriangle(ref direction);
            }

            return true;
        }

        private bool ContainsOriginTriangle(ref Vector3 direction)
        {
            //порядок точек такой, что a всегда последняя добавленная точка
            SupportPoint a = points[0];
            SupportPoint b = points[1];
            SupportPoint c = points[2];

            Vector3 ab = b.Point - a.Point;
            Vector3 ac = c.Point - a.Point;
            Vector3 ao = -a.Point;

            Vector3 abc = Vector3.Cross(ab, ac);

            if (Vector3.Dot(Vector3.Cross(abc, ac), ao) > 0)
            {
                if (Vector3.Dot(ac, ao) > 0)
                {
                    points.Clear();
                    points.Add(a);
                    points.Add(c);
                    //points = { a, c };
                    direction = Vector3.Cross(Vector3.Cross(ac, ao), ac);
                }
                else
                {
                    points.Clear();
                    points.Add(a);
                    points.Add(b);
                    //points = { a, b };
                    return ContainsOriginLine(ref direction);
                }
            }

            else
            {
                if (Vector3.Dot(Vector3.Cross(ab, abc), ao) > 0)
                {
                    points.Clear();
                    points.Add(a);
                    points.Add(b);
                    //points = { a, b }
                    return ContainsOriginLine(ref direction);
                }

                else
                {
                    if (Vector3.Dot(abc, ao) > 0)
                    {
                        direction = abc;
                    }

                    else
                    {
                        points.Clear();
                        points.Add(a);
                        points.Add(c);
                        points.Add(b);
                        //points = { a, c, b };
                        direction = -abc;
                    }
                }
            }

            return false;
        }

        private bool ContainsOriginLine(ref Vector3 direction)
        {
            //порядок точек такой, что a всегда вновь добавленная точка
            SupportPoint a = points[0];
            SupportPoint b = points[1];

            Vector3 ab = b.Point - a.Point;

            Vector3 ao = -a.Point;

            if (Vector3.Dot(ab, ao) > 0)
                //строим перпендикуляр к линии в направлении начала координат
                direction = Vector3.Cross(Vector3.Cross(ab, ao), ab);
            else
            {
                points.Clear();
                points.Add(a);

                direction = ao;
            }

            return false;
        }

    }

#if GJKEPA_DEBUG                
    private static void DrawPolytope(List<SupportPoint> poly, List<int> faces)
    {
        for (int i = 0; i < faces.Count; i += 3)
        {
            Vector3 a = poly[faces[i]].Point;
            Vector3 b = poly[faces[i + 1]].Point;
            Vector3 c = poly[faces[i + 2]].Point;

            UInt32 clr = 4294967040u; //Yellow
            DrawLine(a, b, clr);
            DrawLine(a, c, clr);
            DrawLine(c, b, clr);
        }
    }
#endif

    private static (List<Vector4>, int) GetFaceNormals(List<SupportPoint> polytope, List<int> faces)
    {
        List<Vector4> normals = new List<Vector4>();
        int minTriangle = 0;
        float minDistance = float.MaxValue;

        for (int i = 0; i < faces.Count; i += 3)
        {
            Vector3 a = polytope[faces[i]].Point;
            Vector3 b = polytope[faces[i + 1]].Point;
            Vector3 c = polytope[faces[i + 2]].Point;

            Vector3 normal = Vector3.Cross(b - a, c - a);
            float l = normal.Length();
            float distance = float.MaxValue;

            //If vectors is colliniar we have degenerate face
            if (l < Epsilon)
            {
                normal = Vector3.Zero;
                distance = float.MaxValue;
            }
            else
            {
                normal = normal / l;
                distance = Vector3.Dot(normal, a);
            }

            //If origin outside the polytope
            if (distance < 0)
            {
                normal *= -1;
                distance *= -1;
            }

            normals.Add(new Vector4(normal, distance));

            if (distance < minDistance)
            {
                minTriangle = i / 3;
                minDistance = distance;
            }
        }

        return (normals, minTriangle);
    }

    private static void AddIfUniqueEdge(ref List<(int, int)> edges, List<int> faces, int a, int b)
    {
        //      0--<--3
        //     / \ B /   A: 2-0
        //    / A \ /    B: 0-2
        //   1-->--2
        (int, int) reverse = (0, 0);
        bool found = false;
        foreach ((int, int) f in edges)
        {
            if ((f.Item1 == faces[b]) && (f.Item2 == faces[a]))
            {
                reverse = f;
                found = true;
                break;
            }
        }

        if (found)
        {
            edges.Remove(reverse);
        }
        else
        {
            edges.Add((faces[a], faces[b]));
        }
    }

    private static Vector3 V4to3(Vector4 v)
    {
        return new Vector3(v.X, v.Y, v.Z);
    }

    private static int EPA(Vector3[] polyhedron1, Vector3[] polyhedron2, Simplex simplex, out Vector3 contactNormal, out float PenetrationDepth, out Vector3 contactPoint)
    {
        int epa_cnt = 0;

        List<SupportPoint> polytope = simplex.Vertices;

        //3 vertex indexes for each face
        List<int> faces = new List<int>{
            0, 1, 2,
            0, 3, 1,
            0, 2, 3,
            1, 3, 2
        };

        //we do not know the signs of the projections of the origin onto the normals, but they must be the same
        var (normals, minFace) = GetFaceNormals(polytope, faces);

        Vector3 minNormal = V4to3(normals[minFace]);
        float minDistance = float.MaxValue;

        while ((minDistance == float.MaxValue) && (epa_cnt < epa_max))
        {
            minNormal = V4to3(normals[minFace]);
            minDistance = normals[minFace].W;

            SupportPoint support = CalculateSupportPoint(polyhedron1, polyhedron2, minNormal);
            float sDistance = Vector3.Dot(minNormal, support.Point);

            if (Math.Abs(sDistance - minDistance) > Epsilon)
            {
                minDistance = float.MaxValue;

                List<(int, int)> uniqueEdges = new List<(int, int)>();

                for (int i = 0; i < normals.Count; i++)
                {
                    if ((Vector3.Dot(V4to3(normals[i]), support.Point) - normals[i].W) > 0)
                    {
                        int f = i * 3;
                        AddIfUniqueEdge(ref uniqueEdges, faces, f, f + 1);
                        AddIfUniqueEdge(ref uniqueEdges, faces, f + 1, f + 2);
                        AddIfUniqueEdge(ref uniqueEdges, faces, f + 2, f);
                        faces[f + 2] = faces.Last<int>(); faces.RemoveAt(faces.Count - 1);
                        faces[f + 1] = faces.Last<int>(); faces.RemoveAt(faces.Count - 1);
                        faces[f] = faces.Last<int>(); faces.RemoveAt(faces.Count - 1);

                        normals[i] = normals.Last<Vector4>(); normals.RemoveAt(normals.Count - 1);

                        i--;
                    }
                }

                //!Возможна ситуация при которой из-за невалидных точек (нули, NaN) uniqueEdges не содержит ребер, надо проработать                

                List<int> newFaces = new List<int>();
                for (int i = 0; i < uniqueEdges.Count; i++)
                {
                    (int, int) f = uniqueEdges[i];
                    newFaces.Add(f.Item1);
                    newFaces.Add(f.Item2);
                    newFaces.Add(polytope.Count);

                    //направление соблюдается против часовой
                }

                polytope.Add(support);

                var (newNormals, newMinFace) = GetFaceNormals(polytope, newFaces);

                float oldMinDistance = float.MaxValue;
                for (int i = 0; i < normals.Count; i++)
                {
                    if (normals[i].W < oldMinDistance)
                    {
                        oldMinDistance = normals[i].W;
                        minFace = i;
                    }
                }

                if (newNormals[newMinFace].W < oldMinDistance)
                {
                    minFace = newMinFace + normals.Count;
                }

                faces.AddRange(newFaces);
                normals.AddRange(newNormals);
            }
            epa_cnt++;
        }


#if GJKEPA_DEBUG   
        /*
        DrawPolytope(polytope, faces);
        UInt32 clr = 4294901760u; // Red
        Vector3[] minTri = new Vector3[3];
        minTri[0] = polytope[faces[minFace * 3]].Point;
        minTri[1] = polytope[faces[minFace * 3 + 1]].Point;
        minTri[2] = polytope[faces[minFace * 3 + 2]].Point;
        DrawLine(minTri[0] * 1.1f, minTri[1] * 1.1f, clr);
        DrawLine(minTri[1] * 1.1f, minTri[2] * 1.1f, clr);
        DrawLine(minTri[0] * 1.1f, minTri[2] * 1.1f, clr);
        */
#endif

        SupportPoint a = polytope[faces[minFace * 3]];
        SupportPoint b = polytope[faces[minFace * 3 + 1]];
        SupportPoint c = polytope[faces[minFace * 3 + 2]];

        //Находим проекцию начала координат на плоскость треугольника
        float distance = Vector3.Dot(a.Point, minNormal);
        Vector3 projectedPoint = -distance * minNormal;

        //Получаем барицентрические координаты этой проекции в треугольнике, принадлежащем симпексу
        (float u, float v, float w) = GetBarycentricCoordinates(projectedPoint, a.Point, b.Point, c.Point);

        //Берем соотвествующий треугольник из 1 многогранника
        Vector3 a1 = polyhedron1[a.Index1];
        Vector3 b1 = polyhedron1[b.Index1];
        Vector3 c1 = polyhedron1[c.Index1];

        //Точка контакта на 1 многограннике
        Vector3 contactPoint1 = u * a1 + v * b1 + w * c1;

        //Берем соотвествующий треугольник из 2 многогранника
        Vector3 a2 = polyhedron2[a.Index2];
        Vector3 b2 = polyhedron2[b.Index2];
        Vector3 c2 = polyhedron2[c.Index2];

        //Точка контакта на 2 многограннике
        Vector3 contactPoint2 = u * a2 + v * b2 + w * c2;

        contactPoint = (contactPoint1 + contactPoint2) / 2; //Возвращаем среднюю точку
        contactNormal = minNormal;
        PenetrationDepth = minDistance + 0.001f;

        return epa_cnt;
    }

    public static (float u, float v, float w) GetBarycentricCoordinates(Vector3 p, Vector3 a, Vector3 b, Vector3 c)
    {
        // Векторы от вершины A к вершинам B и C
        Vector3 v0 = b - a, v1 = c - a, v2 = p - a;

        // Вычисляем скалярные произведения
        float d00 = Vector3.Dot(v0, v0);
        float d01 = Vector3.Dot(v0, v1);
        float d11 = Vector3.Dot(v1, v1);
        float d20 = Vector3.Dot(v2, v0);
        float d21 = Vector3.Dot(v2, v1);
        float denom = d00 * d11 - d01 * d01;

        // Check for a zero denominator before division
        if (Math.Abs(denom) <= Epsilon)
        {
            throw new InvalidOperationException("Cannot compute barycentric coordinates for a degenerate triangle.");
        }

        // Вычисляем барицентрические координаты
        float v = (d11 * d20 - d01 * d21) / denom;
        float w = (d00 * d21 - d01 * d20) / denom;
        float u = 1.0f - v - w;

        return (u, v, w);
    }
  
}
