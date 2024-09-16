using System;
using System.Collections;
using System.Collections.Generic;
using System.Windows.Forms;

using Tekla.Structures.Geometry3d;
using Tekla.Structures.Model;
using Tekla.Structures.Model.UI;
using Tekla.Structures.Plugins;
using static Tekla.Structures.Model.Position;
using Tekla.Structures.Solid;
using Identifier = Tekla.Structures.Identifier;

namespace TeklaModelPlugin1

{
    public class PluginData
    {
        #region Fields
        [StructuresField("Thickness")]
        public double Thickness;
        [StructuresField("Material")]
        public string Material;
        [StructuresField("Width")]
        public double Width;
        [StructuresField("DistenceBtwPlates")]
        public double DistenceBtwPlates;
        [StructuresField("NumberOfPlates")]
        public int NumberOfPlates;
        [StructuresField("Name")]
        public string Name;
        [StructuresField("Finish")]
        public string Finish;
        [StructuresField("Offset1")]
        public double Offset1;
        [StructuresField("Offset2")]
        public double Offset2;
        [StructuresField("Reverse")]
        public int Reverse;
        [StructuresField("LayoutIndex")]
        public int LayoutIndex;
        [StructuresField("PlateDepth")]
        public int PlateDepth;
        [StructuresField("Position")]
        public int Position;
        [StructuresField("offsetPosition")]
        public int offsetPosition;
        #endregion
    }

    [Plugin("Batten_Plate_1221")]
    [PluginUserInterface("TeklaModelPlugin1.MainForm")]
    
    public class TeklaModelPlugin1 : PluginBase
    {
        Model myModel = new Model();

        #region Fields
        private Model _Model;
        private PluginData _Data;
        public double _Thickness;
        public string _Material;
        public double _Width;
        public double _DistenceBtwPlates;
        public int _NumberOfPlates;
        public string _Name;
        public string _Finish;
        public double _Offset1;
        public double _Offset2;
        public int _Reverse;
        public int _LayoutIndex;
        public int _PlateDepth;
        public int _Position;
        public int _offsetPosition;
        #endregion

        #region Properties
        private Model Model
        {
            get { return this._Model; }
            set { this._Model = value; }
        }

        private PluginData Data
        {
            get { return this._Data; }
            set { this._Data = value; }
        }
        #endregion

        #region Constructor
        public TeklaModelPlugin1(PluginData data)
        {
            Model = new Model();
            Data = data;
        }
        #endregion

        #region Overrides
        public override List<InputDefinition> DefineInput()
        {
            List<InputDefinition> input = new List<InputDefinition>();
            try
            { //
              // This is an example for selecting two points; change this to suit your needs.
              //
                
                Picker picker = new Picker();
                var part = picker.PickObject(Picker.PickObjectEnum.PICK_ONE_OBJECT, "Pick one object");
                var partno1 = part;
                input.Add(new InputDefinition(partno1.Identifier));
                Point Point = picker.PickPoint();
                input.Add(new InputDefinition(Point));
                part = picker.PickObject(Picker.PickObjectEnum.PICK_ONE_OBJECT, "Pick one object");
                var partno2 = part;
                input.Add(new InputDefinition(partno2.Identifier));

                return input;
            }
            catch { }
            return input;
        }



        public override bool Run(List<InputDefinition> Input)
        {
            try
            {
                GetValuesFromDialog();

                try
                {


                    List<Point> bolt_points = new List<Point>();

                    var beam1Input = Input[0];
                    Beam beam1 = myModel.SelectModelObject(beam1Input.GetInput() as Identifier) as Beam;
                    Point point = Input[1].GetInput() as Point;
                    Beam beam2 = myModel.SelectModelObject(Input[2].GetInput() as Identifier) as Beam;
                    ArrayList beam1_centerpoints = beam1.GetCenterLine(false);
                    ArrayList beam2_centerpoints = beam2.GetCenterLine(false);
                    Point start_point = GetClosestPointOnLineSegment(point, beam1_centerpoints[0] as Point, beam1_centerpoints[1] as Point);
                    Point refference = (_Reverse == 0) ? beam1_centerpoints[1] as Point : beam1_centerpoints[0] as Point;
                    List<Face_> beam1_faces = get_faces(beam1);
                    List<Face_> beam2_faces = get_faces(beam2);
                    Point p1 = new Point(), p2 = new Point(), p3 = new Point(), p4 = new Point(), hold = start_point;
                    List<ContourPlate> plates = new List<ContourPlate>();
                    Line line1 = new Line() , line2 = new Line();
                    switch (_offsetPosition)
                    {
                        case 0:
                            line1 = null;line2 = null;break;
                        case 1:
                            line1 = Intersection.PlaneToPlane(ConvertFaceToGeometricPlane(beam1_faces[(_Position == 1)? 4:0].Face), ConvertFaceToGeometricPlane(beam1_faces[(_Position == 1) ? 3 : 1].Face));
                            line2 = Intersection.PlaneToPlane(ConvertFaceToGeometricPlane(beam1_faces[(_Position == 1) ? 3 : 1].Face), ConvertFaceToGeometricPlane(beam2_faces[(_Position == 1) ? 4 : 0].Face));
                            break;
                        case 2:
                            line1 = Intersection.PlaneToPlane(ConvertFaceToGeometricPlane(beam1_faces[ 2].Face), ConvertFaceToGeometricPlane(beam1_faces[(_Position == 1) ? 3 : 1].Face));
                            line2 = Intersection.PlaneToPlane(ConvertFaceToGeometricPlane(beam1_faces[(_Position == 1) ? 3 : 1].Face), ConvertFaceToGeometricPlane(beam2_faces[2].Face));
                            break ;

                    }


                    if (_LayoutIndex == 0 || _LayoutIndex == 1)
                    {
                        int positionIndex = (_Position == 0) ? 1 : 3;
                        GeometricPlane plane = ConvertFaceToGeometricPlane(beam1_faces[positionIndex].Face);
                        hold = FindPointOnLine(start_point, refference, _DistenceBtwPlates);
                        for (int i = 0; i < _NumberOfPlates; i++)
                        {
                            if (Distance.PointToPoint(hold, refference) >= _Width)//doesnot alows the plates to be created out side the beam
                            {
                                p1 = hold;
                                p2 = FindPointOnLine(p1, refference, _Width);
                                p3 = FindPerpendicularIntersection(beam1_centerpoints[0] as Point, beam1_centerpoints[1] as Point, p2, beam2_centerpoints[0] as Point, beam2_centerpoints[1] as Point);
                                p4 = FindPerpendicularIntersection(beam1_centerpoints[0] as Point, beam1_centerpoints[1] as Point, p1, beam2_centerpoints[0] as Point, beam2_centerpoints[1] as Point);

                                plates.Add(countourPlate(FindClosestPointOnPlane(plane, p1), FindClosestPointOnPlane(plane, p2), FindClosestPointOnPlane(plane, p3), FindClosestPointOnPlane(plane, p4), beam1, beam2, false, line1, line2));
                                hold = FindPointOnLine(p2, refference, _DistenceBtwPlates);
                            }
                        }
                    }
                    hold = start_point;
                    switch (_offsetPosition)
                    {
                        case 0:
                            line1 = null; line2 = null; break;
                        case 1:
                            line1 = Intersection.PlaneToPlane(ConvertFaceToGeometricPlane(beam1_faces[(_Position == 1) ? 4 : 0].Face), ConvertFaceToGeometricPlane(beam1_faces[(_Position == 1) ? 5 : 7].Face));
                            line2 = Intersection.PlaneToPlane(ConvertFaceToGeometricPlane(beam1_faces[(_Position == 1) ? 5 : 7].Face), ConvertFaceToGeometricPlane(beam2_faces[(_Position == 1) ? 4 : 0].Face));
                            break;
                        case 2:
                            line1 = Intersection.PlaneToPlane(ConvertFaceToGeometricPlane(beam1_faces[6].Face), ConvertFaceToGeometricPlane(beam1_faces[(_Position == 1) ? 5 : 7].Face));
                            line2 = Intersection.PlaneToPlane(ConvertFaceToGeometricPlane(beam1_faces[(_Position == 1) ? 5 : 7].Face), ConvertFaceToGeometricPlane(beam2_faces[2].Face));
                            break;

                    }
                    if (_LayoutIndex == 0 || _LayoutIndex == 2)
                    {
                        hold = FindPointOnLine(start_point, refference, _DistenceBtwPlates);
                        int positionIndex = (_Position == 0) ? 7 : 5;
                        GeometricPlane plane = ConvertFaceToGeometricPlane(beam1_faces[positionIndex].Face);
                        for (int i = 0; i < _NumberOfPlates; i++)
                        {
                            if (Distance.PointToPoint(hold, refference) >= _Width)//doesnot alows the plates to be created out side the beam
                            {
                                p1 = hold;
                                p2 = FindPointOnLine(p1, refference, _Width);
                                p3 = FindPerpendicularIntersection(beam1_centerpoints[0] as Point, beam1_centerpoints[1] as Point, p2, beam2_centerpoints[0] as Point, beam2_centerpoints[1] as Point);
                                p4 = FindPerpendicularIntersection(beam1_centerpoints[0] as Point, beam1_centerpoints[1] as Point, p1, beam2_centerpoints[0] as Point, beam2_centerpoints[1] as Point);
                                plates.Add(countourPlate(FindClosestPointOnPlane(plane, p1), FindClosestPointOnPlane(plane, p2), FindClosestPointOnPlane(plane, p3), FindClosestPointOnPlane(plane, p4), beam1, beam2, true, line1, line2));
                                hold = FindPointOnLine(p2, refference, _DistenceBtwPlates);
                            }
                        }

                    }
                    weld(beam1, beam2, plates);


                    return true;

                }
                catch { }
                return false;
            }
            catch (Exception Exc)
            {
                MessageBox.Show(Exc.ToString());
            }

            return true;
        }
        #endregion

        #region Private methods
        /// <summary>
        /// Gets the values from the dialog and sets the default values if needed
        /// </summary>
        private void GetValuesFromDialog()
        {

            {
                _Thickness = Data.Thickness;
                _Width = Data.Width;
                _Material = Data.Material;
                _DistenceBtwPlates = Data.DistenceBtwPlates;
                _NumberOfPlates = Data.NumberOfPlates;
                _Name = Data.Name;
                _Finish = Data.Finish;
                _Offset1 = Data.Offset1;
                _Offset2 = Data.Offset2;
                _Reverse = Data.Reverse;
                _LayoutIndex = Data.LayoutIndex;
                _PlateDepth = Data.PlateDepth;
                _Position = Data.Position;
                _offsetPosition = Data.offsetPosition;
                if (IsDefaultValue(_Thickness))
                {
                    _Thickness = 10;
                }
                if (IsDefaultValue(_Width))
                {
                    _Width = 300;
                }
                if (IsDefaultValue(_Material))
                {
                    _Material = "IS2062";
                }
                if (IsDefaultValue(_DistenceBtwPlates))
                {
                    _DistenceBtwPlates = 300;
                }
                if (IsDefaultValue(_NumberOfPlates))
                {
                    _NumberOfPlates = 5;
                }
                if (IsDefaultValue(_Name))
                {
                    _Name = "Batten plate";
                }
                if (IsDefaultValue(_Finish))
                {
                    _Finish = "foo";
                }
                if (IsDefaultValue(_Offset1))
                {
                    _Offset1 = 0;
                }
                if (IsDefaultValue(_Offset2))
                {
                    _Offset2 = 0;
                }
                if (IsDefaultValue(_Reverse))
                {
                    _Reverse = 0;

                }
                if (IsDefaultValue(_LayoutIndex))
                {
                    _LayoutIndex = 1;
                }
                if (IsDefaultValue(_PlateDepth))
                {
                    _PlateDepth = 2;
                }
                if (IsDefaultValue(_offsetPosition))
                {
                    _offsetPosition = 0;
                }
            }
        }

        class Face_
        {
            public Face Face { get; set; }
            public Vector Vector { get; set; }
            public void face_(Face face, Vector vector)
            {
                Face = face;
                Vector = vector;
            }
        }
        private List<Face_> get_faces(Beam beam)
        {

            Solid solid = beam.GetSolid();
            FaceEnumerator faceEnumerator = solid.GetFaceEnumerator();
            List<Face_> faces = new List<Face_>();
            while (faceEnumerator.MoveNext())
            {

                Face face = faceEnumerator.Current as Face;
                Vector vector = face.Normal;
                faces.Add(new Face_ { Face = face, Vector = vector });

            }

            return faces;
        }
        private Point midPoint(Point point, Point point1)
        {
            Point mid = new Point((point.X + point1.X) / 2, (point.Y + point1.Y) / 2, (point.Z + point1.Z) / 2);
            return mid;
        }
        public static GeometricPlane ConvertFaceToGeometricPlane(Face face)
        {
            ArrayList points = new ArrayList();
            // Get the edges from the face (since 'Points' is not available)
            LoopEnumerator loopEnumerator = face.GetLoopEnumerator();
            while (loopEnumerator.MoveNext())
            {

                Loop loop = loopEnumerator.Current as Loop;
                VertexEnumerator vertexEnumerator = loop.GetVertexEnumerator();
                while (vertexEnumerator.MoveNext())
                {
                    points.Add(vertexEnumerator.Current);
                }
            }

            Point point1 = points[0] as Point;
            Point point2 = points[1] as Point;
            Point point3 = points[2] as Point;



            if (point1 == null || point2 == null || point3 == null)
            {
                throw new ArgumentException("The face does not have sufficient points to define a plane.");
            }

            // Create vectors from the points
            Vector vector1 = new Vector(point2.X - point1.X, point2.Y - point1.Y, point2.Z - point1.Z);
            Vector vector2 = new Vector(point3.X - point1.X, point3.Y - point1.Y, point3.Z - point1.Z);

            // Calculate the normal vector (cross product of the two vectors)
            Vector normalVector = Vector.Cross(vector1, vector2);
            normalVector.Normalize();

            // Create the geometric plane using point1 and the normal vector
            GeometricPlane geometricPlane = new GeometricPlane(point1, normalVector);

            return geometricPlane;
        }

        private ArrayList get_points(Face face)
        {
            ArrayList points = new ArrayList();
            LoopEnumerator loopEnumerator = face.GetLoopEnumerator();
            while (loopEnumerator.MoveNext())
            {

                Loop loop = loopEnumerator.Current as Loop;
                VertexEnumerator vertexEnumerator = loop.GetVertexEnumerator();
                while (vertexEnumerator.MoveNext())
                {
                    points.Add(vertexEnumerator.Current);
                }
            }
            return points;
        }
        private static double DistanceBetweenPoints(Point point1, Point point2)
        {
            return Math.Sqrt(Math.Pow(point2.X - point1.X, 2) +
                             Math.Pow(point2.Y - point1.Y, 2) +
                             Math.Pow(point2.Z - point1.Z, 2));
        }
        public static Point FindPointOnLine(Point startPoint, Point secondPoint, double distance)
        {
            if (distance == 0)
                return startPoint;
            // Step 1: Calculate the direction vector from startPoint to secondPoint
            Vector direction = new Vector(
                secondPoint.X - startPoint.X,
                secondPoint.Y - startPoint.Y,
                secondPoint.Z - startPoint.Z
            );

            // Step 2: Normalize the direction vector
            direction.Normalize();

            // Step 3: Scale the direction vector by the distance
            Vector scaledVector = new Vector(
                direction.X * distance,
                direction.Y * distance,
                direction.Z * distance
            );

            // Step 4: Calculate the new point by adding the scaled vector to the start point
            Point newPoint = new Point(
                startPoint.X + scaledVector.X,
                startPoint.Y + scaledVector.Y,
                startPoint.Z + scaledVector.Z
            );

            return newPoint;
        }
        private double CalculateDistanceBetweenFaces(Face face1, Face face2)
        {
            // Get the loop vertices of both faces to extract points
            ArrayList face1Vertices = get_points(face1);
            ArrayList face2Vertices = get_points(face2);

            if (face1Vertices == null || face1Vertices.Count == 0 || face2Vertices == null || face2Vertices.Count == 0)
            {
                throw new ArgumentException("One or both faces do not have vertices.");
            }

            // Initialize the minimum distance to a large value
            double minDistance = double.MaxValue;

            // Loop through all points on face1 and face2 and calculate the distance between each pair
            foreach (Point p1 in face1Vertices)
            {
                foreach (Point p2 in face2Vertices)
                {
                    double distance = DistanceBetweenPoints(p1, p2);
                    if (distance < minDistance)
                    {
                        minDistance = distance;
                    }
                }
            }

            return minDistance;
        }
        private static Point GetClosestPointOnLineSegment(Point point, Point lineStart, Point lineEnd)
        {
            // Vector from line start to the point
            Vector startToPoint = new Vector(point - lineStart);

            // Direction vector of the line segment
            Vector lineDirection = new Vector(lineEnd - lineStart);
            double lineLengthSquared = lineDirection.Dot(lineDirection);

            // Project the point onto the line segment
            double t = startToPoint.Dot(lineDirection) / lineLengthSquared;

            // Clamp t to the range [0, 1] to keep the projection within the segment
            t = Math.Max(0, Math.Min(1, t));

            // Calculate the closest point on the line segment
            return new Point(
                lineStart.X + t * lineDirection.X,
                lineStart.Y + t * lineDirection.Y,
                lineStart.Z + t * lineDirection.Z
            );
        }
        private ContourPlate countourPlate(Point p1, Point p2, Point p3, Point p4, Beam beam1, Beam beam2, bool flag ,Line line1,Line line2)
        {
            if (_Reverse == 1)
                flag = !flag;
            if (_Position == 1)
                flag = !flag;
            if (line1 != null && line2 != null)//sets for the edge position
            {
                p1 = Intersection.LineToLine(line1, new Line(p1, p4)).StartPoint;
                p2 = Intersection.LineToLine(line1, new Line(p2, p3)).StartPoint;
                p3 = Intersection.LineToLine(line2, new Line(p3, p2)).StartPoint;
                p4 = Intersection.LineToLine(line2, new Line(p1, p4)).StartPoint;
            }
            p1 = FindPointOnLine(p1, p4, _Offset1);//sets for the offset
            p2 = FindPointOnLine(p2, p3, _Offset1);
            p3 = FindPointOnLine(p3, p2, _Offset2);
            p4 = FindPointOnLine(p4, p1, _Offset2);
           

            ArrayList countp = new ArrayList();

            ContourPoint contourPoint = new ContourPoint(p1, new Chamfer());
            countp.Add(contourPoint);
            contourPoint = new ContourPoint(p2, new Chamfer());
            countp.Add(contourPoint);
            contourPoint = new ContourPoint(p3, new Chamfer());
            countp.Add(contourPoint);
            contourPoint = new ContourPoint(p4, new Chamfer());
            countp.Add(contourPoint);

            ContourPlate cp = new ContourPlate();
            cp.Contour.ContourPoints = countp;

            cp.Profile.ProfileString = "PLT" + _Thickness;

            cp.Material.MaterialString = _Material;
            cp.Class = "4";
            if (_PlateDepth == 2)
                cp.Position.Depth = (flag) ? Position.DepthEnum.FRONT : DepthEnum.BEHIND;
            else
            {
                if (_PlateDepth == 0)
                    cp.Position.Depth = (flag) ? Position.DepthEnum.FRONT : DepthEnum.BEHIND;
                else
                    cp.Position.Depth = (!flag) ? Position.DepthEnum.FRONT : DepthEnum.BEHIND;
            }

            cp.Name = _Name;
            cp.Finish = _Finish;
            cp.Position.DepthOffset = 0;
            cp.Insert();
            

            return cp;

        }

        public static Point FindPerpendicularIntersection(Point line1Start, Point line1End, Point pointOnLine1, Point line2Start, Point line2End)
        {
            Vector vector = new Vector(line1Start.X - line1End.X, line1Start.Y - line1End.Y, line1Start.Z - line1End.Z);
            GeometricPlane plane = new GeometricPlane(pointOnLine1, vector);
            Point intersection = Intersection.LineToPlane(new Line(line2Start, line2End), plane);
            return intersection;
        }

        public static double IsNumber(string input, double default_value)
        {
            // Try to parse as an integer
            if (int.TryParse(input, out int result))
            {
                return result; // It's a valid integer
            }

            // Try to parse as a double
            if (double.TryParse(input, out double result1))
            {
                return result1; // It's a valid double
            }

            // If neither parsing works, it's not a valid number
            return default_value;
        }
        private void weld(Part part1, Part part2, List<ContourPlate> plates)
        {
            try
            {
                Weld Weld = new Weld();
                foreach (ContourPlate p in plates)
                {

                    Weld.MainObject = p;
                    Weld.SecondaryObject = part1;
                    Weld.TypeAbove = BaseWeld.WeldTypeEnum.WELD_TYPE_FILLET;
                    Weld.TypeBelow = BaseWeld.WeldTypeEnum.WELD_TYPE_FILLET;
                    Weld.LengthAbove = 12;
                    Weld.TypeBelow = BaseWeld.WeldTypeEnum.WELD_TYPE_SLOT;
                    Weld.Insert();

                    Weld.SecondaryObject = part2;
                    Weld.TypeAbove = BaseWeld.WeldTypeEnum.WELD_TYPE_FILLET;
                    Weld.TypeBelow = BaseWeld.WeldTypeEnum.WELD_TYPE_FILLET;
                    Weld.LengthAbove = 12;
                    Weld.TypeBelow = BaseWeld.WeldTypeEnum.WELD_TYPE_SLOT;
                    Weld.Insert();
                }

                Weld.Modify();
            }
            catch { }
        }
        public static Point FindClosestPointOnPlane(GeometricPlane plane, Point point)
        {
            // Step 1: Get the normal vector of the plane
            Vector normalVector = plane.Normal;

            // Step 2: Find a vector from the plane's origin to the given point
            Vector pointToPlaneVector = new Vector(point.X - plane.Origin.X, point.Y - plane.Origin.Y, point.Z - plane.Origin.Z);

            // Step 3: Project the pointToPlaneVector onto the plane's normal vector
            double distanceFromPointToPlane = pointToPlaneVector.Dot(normalVector); // Dot product to find projection length along the normal

            // Step 4: Calculate the closest point by moving from the point in the opposite direction of the normal by the distance
            Point closestPoint = new Point(
                point.X - distanceFromPointToPlane * normalVector.X,
                point.Y - distanceFromPointToPlane * normalVector.Y,
                point.Z - distanceFromPointToPlane * normalVector.Z
            );

            return closestPoint;
        }

        #endregion
    }
}
