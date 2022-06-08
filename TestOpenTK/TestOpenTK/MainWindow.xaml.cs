using Common;
using OpenTK.Wpf;
using System;
using System.Runtime.InteropServices;
using System.Windows;
using OpenTK.Graphics.OpenGL4;
using System.Collections.Generic;
using OpenTK;
using System.Windows.Input;
using System.Drawing;
using System.Drawing.Imaging;
using System.Timers;
using System.Windows.Controls;
using System.Windows.Media;
using System.Diagnostics;

namespace TestOpenTK
{
    public class RSFrame
    {
        public float[] PointsX { get; set; }
        public float[] PointsY { get; set; }
        public float[] PointsZ { get; set; }
        public float[] TextureU { get; set; }
        public float[] TextureV { get; set; }
    }

    public class Point3D
    {
        // the original x y z that came from RS before model and view transformations
        public float X { get; set; }
        public float Y { get; set; }
        public float Z { get; set; }
    }

    public class Point3DXY : Point3D
    {
        public float SX { get; set; }
        public float SY { get; set; }
    }

    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        double OpenTKControlWidth = 800;
        double OpenTKControlHeight = 600;
        Timer _timer;
        bool _grabImage = false;
        bool _grabDepthImage = false;
        int _nextImageIndex = 1;
        int _nextDepthIndex = 1;
        bool _cutOutliers = true;

        float _modelRotationX = 0;
        float _modelRotationY = 0;

        float _modelRotationCenterX = 0;
        float _modelRotationCenterY = 0;
        float _modelRotationCenterZ = 0;

        bool _updateRotCenter = true;

        RSFrame _lastFrame; // the last frame arrived from the server
        Dictionary<int, Dictionary<int, List<Point3DXY>>> _xyIndex = null;

        //float _clipx = -1;
        //float _clipy = -1;

        const int XDetla = 3;
        const int YDetla = 3;
        private bool _measureDistance = false;
        //private System.Windows.Point _measureStart;
        private Point3D _measure3DStart = null;


        // This class is a wrapper around a shader, which helps us manage it.
        // The shader class's code is in the Common project.
        // What shaders are and what they're used for will be explained later in this tutorial.
        private Shader _shader;
        private Bitmap _bmp = null;

        // We need an instance of the new camera class so it can manage the view and projection matrix code.
        // We also need a boolean set to true to detect whether or not the mouse has been moved for the first time.
        // Finally, we add the last position of the mouse so we can calculate the mouse offset easily.
        private Camera _camera;

        private int _vertexBufferObject;

        private int _vertexArrayObject;

        // Create the vertices for our triangle. These are listed in normalized device coordinates (NDC)
        // In NDC, (0, 0) is the center of the screen.
        // Negative X coordinates move to the left, positive X move to the right.
        // Negative Y coordinates move to the bottom, positive Y move to the top.
        // OpenGL only supports rendering in 3D, so to create a flat triangle, the Z coordinate will be kept as 0.
        private float[] _vertices = { };
        private object _syncLock = new object();

        private Texture _texture;

        bool _isCapturing = false;

        [DllImport("rs-pointcloud.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern bool GetNextFrame(
            out int textureFormat,
            out int textureWidth,
            out int textureHeight,
            out IntPtr texturePixels,
            out int pointsSize,
            out IntPtr pointsX,
            out IntPtr pointsY,
            out IntPtr pointsZ,
            out IntPtr textureU,
            out IntPtr textureV);

        [DllImport("rs-pointcloud.dll", CallingConvention = CallingConvention.Cdecl)]
        public static extern void Initialize();


        bool executing = false;
        private void _timer_Elapsed(object sender, ElapsedEventArgs e)
        {
            if (!_isCapturing || executing)
            {
                return;
            }

            executing = true;
//            Dispatcher.Invoke(new Action(() =>
  //          {
                GrabFrame();
            //}));

            System.Threading.Thread.Sleep(1);
            executing = false;
        }

        public MainWindow()
        {
            InitializeComponent();
            var settings = new GLWpfControlSettings
            {
                MajorVersion = 3,
                MinorVersion = 6
            };

            OpenTkControl.Start(settings);
            Initialize();

            _timer = new Timer(10);
            _timer.Elapsed += _timer_Elapsed;
            _timer.Start();
        }

        private void OpenTkControl_OnRender(TimeSpan delta)
        {
            lock (_syncLock)
            {
                Init();
                GL.Clear(ClearBufferMask.ColorBufferBit);

                var vertices = _vertices;
                var bmp = _bmp;
                _bmp = null;

                GL.BufferData(BufferTarget.ArrayBuffer, vertices.Length * sizeof(float), vertices, BufferUsageHint.StaticDraw);

                if (bmp != null)
                {
                    using (bmp)
                    {
                        _texture = Texture.LoadFromBitmap(bmp);
                        _texture.Use(TextureUnit.Texture0);
                    }
                }

                // Bind the shader
                _shader.Use();

                // Bind the VAO
                GL.BindVertexArray(_vertexArrayObject);

                var model = Matrix4.Identity * Matrix4.CreateRotationX(0);
                _shader.SetMatrix4("model", model);
                _shader.SetMatrix4("view", _camera.GetViewMatrix());
                _shader.SetMatrix4("projection", _camera.GetProjectionMatrix());

                GL.DrawArrays(PrimitiveType.Points, 0, vertices.Length);
            }
        }

        private Matrix3? GetModelRotationMatrix()
        {
            if (_modelRotationX == 0d && _modelRotationY == 0d)
            {
                return null;
            }

            Matrix3 rotx;
            Matrix3.CreateRotationX(_modelRotationX, out rotx);
            Matrix3 roty;
            Matrix3.CreateRotationY(_modelRotationY, out roty);
            return rotx * roty;
        }

        private void ApplyModelRotation(ref float x, ref float y, ref float z, Matrix3 rot)
        {
            var vec = new Vector3(x- _modelRotationCenterX, y - _modelRotationCenterY, z - _modelRotationCenterZ);
            Vector3 res;
            Vector3.Transform(ref vec, ref rot, out res);
            x = res.X + _modelRotationCenterX;
            y = res.Y + _modelRotationCenterY;
            z = res.Z + _modelRotationCenterZ;
        }

        /// <summary>
        /// Updates _xyIndex, _modelRotationCenterX,Y,Z and returns the result Vertex Array
        /// </summary>
        /// <param name="rsFrame"></param>
        /// <returns></returns>
        private float[] BuildVertexArray(RSFrame rsFrame)
        {
            decimal xSum = 0;
            decimal ySum = 0;
            decimal zSum = 0;
            int numVals = 0;
            bool updateRot = _updateRotCenter;

            var modelRotation = GetModelRotationMatrix();
            var xyIndex = new Dictionary<int, Dictionary<int, List<Point3DXY>>>();
            List<float> res = new List<float>();
            for (int i = 0; i < rsFrame.PointsZ.Length; ++i)
            {
                if (rsFrame.PointsZ[i] != 0 && (!_cutOutliers ||
                    (rsFrame.TextureU[i] > 0 && rsFrame.TextureU[i] <= 1
                    && rsFrame.TextureV[i] > 0 && rsFrame.TextureV[i] <= 1)))
                {
                    var originalX = rsFrame.PointsX[i];
                    var originalY = -rsFrame.PointsY[i];
                    var originalZ = -rsFrame.PointsZ[i];

                    var x = originalX;
                    var y = originalY;
                    var z = originalZ;

                    if (updateRot)
                    {
                        xSum += Convert.ToDecimal(x);
                        ySum += Convert.ToDecimal(y);
                        zSum += Convert.ToDecimal(z);
                        ++numVals;
                    }

                    if (modelRotation.HasValue)
                    {
                        ApplyModelRotation(ref x, ref y, ref z, modelRotation.Value);
                    }

                    var w = Convert.ToSingle(OpenTKControlWidth);
                    var h = Convert.ToSingle(OpenTKControlHeight);

                    float sx, sy;
                    GetScreenPos(x, y, z, 
                        //_clipx, _clipy,
                        w, h, out sx, out sy);

                    var xInd = Convert.ToInt32(Math.Floor(sx / XDetla));
                    var yInd = Convert.ToInt32(Math.Floor(sy / YDetla));

                    // TODO: test 0.7/3 in double
                    //const int YDetla = 3;
                    Dictionary<int, List<Point3DXY>> xDict;
                    if (!xyIndex.TryGetValue(xInd, out xDict))
                    {
                        xDict = new Dictionary<int, List<Point3DXY>>();
                        xyIndex[xInd] = xDict;
                    }

                    List<Point3DXY> yList;
                    if (!xDict.TryGetValue(yInd, out yList))
                    {
                        yList = new List<Point3DXY>();
                        xDict[yInd] = yList;
                    }

                    yList.Add(new Point3DXY
                    {
                        X = originalX,
                        Y = originalY,
                        Z = originalZ,
                        SX = sx,
                        SY = sy
                    });

                    //if (d > 10)
                    //{
                    res.Add(x);
                    res.Add(y);
                    res.Add(z);
                    res.Add(rsFrame.TextureU[i]);
                    res.Add(rsFrame.TextureV[i]);
                    //}
                }
            }
            _xyIndex = xyIndex;

            // update measuring line start
            if (_measure3DStart != null)
            {
                var x = _measure3DStart.X;
                var y = _measure3DStart.Y;
                var z = _measure3DStart.Z;
                Debug.WriteLine($"x:{x},y:{y},z:{z},");

                if (modelRotation.HasValue)
                {
                    ApplyModelRotation(ref x, ref y, ref z, modelRotation.Value);
                }

                var w = Convert.ToSingle(OpenTKControlWidth);
                var h = Convert.ToSingle(OpenTKControlHeight);

                float sx, sy;
                GetScreenPos(x, y, z,
                    //_clipx, _clipy,
                    w, h, out sx, out sy);
                Dispatcher.BeginInvoke(new Action(() =>
                {
                    distanceLine.X1 = sx;
                    distanceLine.Y1 = sy;
                }));
            }

            if (updateRot && numVals > 0)
            {
                _modelRotationCenterX = Convert.ToSingle(xSum/ numVals);
                _modelRotationCenterY = Convert.ToSingle(ySum / numVals);
                _modelRotationCenterZ = Convert.ToSingle(zSum / numVals);
            }
            return res.ToArray();
        }

        private Bitmap BMPFromIntPtr(IntPtr intptr, int size, int width, int height)
        {
            byte[] buffer = new byte[size];
            Marshal.Copy(intptr, buffer, 0, size);

            Bitmap b = new Bitmap(width, height, System.Drawing.Imaging.PixelFormat.Format24bppRgb);

            Rectangle BoundsRect = new Rectangle(0, 0, width, height);
            BitmapData bmpData = b.LockBits(BoundsRect,
                                            ImageLockMode.WriteOnly,
                                            b.PixelFormat);

            IntPtr ptr = bmpData.Scan0;

            // add back dummy bytes between lines, make each line be a multiple of 4 bytes
            int skipByte = bmpData.Stride - width * 3;
            byte[] newBuff = new byte[buffer.Length + skipByte * height];
            for (int j = 0; j < height; j++)
            {
                System.Buffer.BlockCopy(buffer, j * width * 3, newBuff, j * (width * 3 + skipByte), width * 3);
            }

            // fill in rgbValues
            Marshal.Copy(newBuff, 0, ptr, newBuff.Length);
            b.UnlockBits(bmpData);
            if (_grabImage)
            {
                var fileName = @"c:\Misc\myPic" + _nextImageIndex.ToString() + ".bmp";
                ++_nextImageIndex;
                if (!System.IO.File.Exists(fileName))
                {
                    b.Save(fileName, ImageFormat.Bmp);
                }
                _grabImage = false;
            }
            return b;
        }

        //private static void MySaveBMP(IntPtr intptr, int size, int width, int height)
        //{
        //    byte[] buffer = new byte[size];
        //    Marshal.Copy(intptr, buffer, 0, size);

        //    Bitmap b = new Bitmap(width, height, System.Drawing.Imaging.PixelFormat.Format24bppRgb);

        //    Rectangle BoundsRect = new Rectangle(0, 0, width, height);
        //    BitmapData bmpData = b.LockBits(BoundsRect,
        //                                    ImageLockMode.WriteOnly,
        //                                    b.PixelFormat);

        //    IntPtr ptr = bmpData.Scan0;

        //    // add back dummy bytes between lines, make each line be a multiple of 4 bytes
        //    int skipByte = bmpData.Stride - width * 3;
        //    byte[] newBuff = new byte[buffer.Length + skipByte * height];
        //    for (int j = 0; j < height; j++)
        //    {
        //        System.Buffer.BlockCopy(buffer, j * width * 3, newBuff, j * (width * 3 + skipByte), width * 3);
        //    }

        //    // fill in rgbValues
        //    Marshal.Copy(newBuff, 0, ptr, newBuff.Length);
        //    b.UnlockBits(bmpData);
        //    b.Save(@"c:\Misc\myPic.bmp", ImageFormat.Bmp);
        //}

        class Depthcand
        {
            public double PosX { get; set; }
            public double PosY { get; set; }
            public float Depth { get; set; }
        }

        //public class ImageCoords
        //{
        //    List<>
        //}

        //private float[] GetDepthImageMatrix(
        //    int width,
        //    int height,
        //    float[] pointsX,
        //    float[] pointsY,
        //    float[] pointsZ,
        //    float[] textureU,
        //    float[] textureV,
        //    out byte[] depthImage)
        //{
        //    var depthCands = new Dictionary<int, List<Depthcand>>();
        //    var minDepth = float.MaxValue;
        //    var maxDepth = float.MinValue;
        //    for (int pci = 0; pci < pointsX.Length; ++pci)
        //    {
        //        var px = pointsX[pci];
        //        var py = pointsY[pci];
        //        var pz = pointsZ[pci];
        //        float distance = Convert.ToSingle(Math.Sqrt(px * px + py * py + pz * pz));
        //        if (distance < minDepth)
        //        {
        //            minDepth = distance;
        //        }
        //        if (distance > maxDepth)
        //        {
        //            maxDepth = distance;
        //        }

        //        var uwidth = textureU[pci] * width;
        //        uwidth = uwidth < 0 ? 0 : uwidth;
        //        uwidth = uwidth > width - 1 ? width - 1 : uwidth;
        //        var vheight = textureV[pci] * height;
        //        vheight = vheight < 0 ? 0 : vheight;
        //        vheight = vheight > height - 1 ? height - 1 : vheight - 1;
        //        int rwidth = Convert.ToInt32(Math.Round(uwidth));
        //        int rheight = Convert.ToInt32(Math.Round(vheight));
        //        List<Depthcand> cands;
        //        var ind = rwidth * rheight + rwidth;
        //        if (!depthCands.TryGetValue(ind, out cands))
        //        {
        //            cands = new List<Depthcand>();
        //            depthCands[ind] = cands;
        //        }

        //        cands.Add(new Depthcand { PosX= textureU[pci] * width, PosY = textureV[pci] * height, Depth = distance });
        //    }

        //    const float undefinedValue = 0;
        //    var result = new float[width * height];
        //    depthImage = new byte[width * height];

        //    int delta = 2;
        //    double maxd = 2d;
        //    double maxd2 = maxd* maxd;
        //    for (int x = 0; x < width; ++x)
        //    {
        //        for (int y = 0;y < height; ++y)
        //        {
        //            Depthcand bestCand = null;
        //            var bestDistance2 = maxd2;
        //            var bestDepth = undefinedValue;
        //            for (int dx = -delta; dx <= delta; ++dx)
        //            {
        //                for (int dy = -delta; dy <= delta; ++dy)
        //                {
        //                    var xpos = x + dx;
        //                    var ypos = x + dx;
        //                    if (xpos < 0 || xpos >= width || ypos < 0 || ypos >= height)
        //                    {
        //                        continue;
        //                    }
        //                    var ind = ypos * width + xpos;
        //                    List<Depthcand> candsList;
        //                    if (depthCands.TryGetValue(ind, out candsList))
        //                    {
        //                        foreach (var cand in candsList)
        //                        {
        //                            var ddx = cand.PosX - x;
        //                            var ddy = cand.PosY - y;
        //                            var d2 = ddx * ddx + ddy * ddy;
        //                            if (d2 < bestDistance2)
        //                            {
        //                                bestCand = cand;
        //                                bestDistance2 = d2;
        //                                bestDepth =  cand.Depth;
        //                            }
        //                        }
        //                    }
        //                }
        //            }

        //            result[y * width + x] = bestDepth;
        //            depthImage[y * width + x] =
        //                (maxDepth == minDepth || bestDepth == undefinedValue)  ? Convert.ToByte(0) :
        //                Convert.ToByte(((bestDepth - minDepth) / (maxDepth - minDepth)) * 200 + 55);
        //        }
        //    }

        //    return result;
        //}

        private void DepthToImage(int width, int height)
        {
            var b = new Bitmap(width, height, System.Drawing.Imaging.PixelFormat.Format8bppIndexed);

            ColorPalette ncp = b.Palette;
            for (int i = 0; i < 256; i++)
            {
                ncp.Entries[i] = System.Drawing.Color.FromArgb(255, i, i, i);
            }
            b.Palette = ncp;

            var BoundsRect = new Rectangle(0, 0, width, height);
            BitmapData bmpData = b.LockBits(BoundsRect,
                                            ImageLockMode.WriteOnly,
                                            b.PixelFormat);

            IntPtr ptr = bmpData.Scan0;

            int bytes = bmpData.Stride * b.Height;
            var rgbValues = new byte[bytes];
            for (int w = 0; w < bmpData.Stride; ++w)
            {
                for (int h = 0; h < b.Height; ++h)
                {
                    rgbValues[h * bmpData.Stride + w] = Convert.ToByte(255 * ((bmpData.Stride - w) * (b.Height - h)) / (bmpData.Stride * b.Height));
                }
            }

            Marshal.Copy(rgbValues, 0, ptr, bytes);
            b.UnlockBits(bmpData);
            if (_grabDepthImage)
            {
                var fileName = @"c:\Misc\myDepthPic" + _nextDepthIndex.ToString() + ".bmp";
                ++_nextDepthIndex;
                if (!System.IO.File.Exists(fileName))
                {
                    b.Save(fileName, ImageFormat.Bmp);
                }
                _grabDepthImage = false;
            }
        }

        private void UpdateNextView(Bitmap bmp,
            RSFrame rSFrame)
        {
            lock (_syncLock)
            {
                _lastFrame = rSFrame;
                _bmp = bmp;
            }
        }

        private void GrabFrame()
        {
            int textureFormatOut;
            int textureWidthOut;
            int textureHeightOut;
            IntPtr texturePixelsOut;
            int pointsSizeOut;
            IntPtr pointsXOut;
            IntPtr pointsYOut;
            IntPtr pointsZOut;
            IntPtr textureUOut;
            IntPtr textureVOut;

            try
            {
                if (!GetNextFrame(
                    out textureFormatOut,
                    out textureWidthOut,
                    out textureHeightOut,
                    out texturePixelsOut,
                    out pointsSizeOut,
                    out pointsXOut,
                    out pointsYOut,
                    out pointsZOut,
                    out textureUOut,
                    out textureVOut))
                {
                    return;
                }
            }
            catch (Exception ex)
            {
                MessageBox.Show("Is Realsense Camera connected? if not, connect the camera and restart the application.");
                this.Close();
                _timer.Stop();
                _timer.Dispose();
                return;
            }

            Bitmap bmp = BMPFromIntPtr(texturePixelsOut, textureWidthOut * textureHeightOut * 3, textureWidthOut, textureHeightOut);
            if (_grabDepthImage) // save the image to a file if requested
            {
                DepthToImage(textureWidthOut, textureHeightOut);
            }

            float[] pointsX = new float[pointsSizeOut];
            float[] pointsY = new float[pointsSizeOut];
            float[] pointsZ = new float[pointsSizeOut];
            float[] textureU = new float[pointsSizeOut];
            float[] textureV = new float[pointsSizeOut];

            Marshal.Copy(pointsXOut, pointsX, 0, pointsSizeOut);
            Marshal.Copy(pointsYOut, pointsY, 0, pointsSizeOut);
            Marshal.Copy(pointsZOut, pointsZ, 0, pointsSizeOut);
            Marshal.Copy(textureUOut, textureU, 0, pointsSizeOut);
            Marshal.Copy(textureVOut, textureV, 0, pointsSizeOut);

            var resFrame = new RSFrame
            {
                PointsX = pointsX,
                PointsY = pointsY,
                PointsZ = pointsZ,
                TextureU = textureU,
                TextureV = textureV
            };

            UpdateNextView(bmp, resFrame);
            _vertices = BuildVertexArray(_lastFrame);
            Dispatcher.Invoke(RefreshView);
        }

        private void RefreshView()
        {
            lock (_syncLock)
            {
                if (_lastFrame == null)
                {
                    return;
                }

                //_vertices = BuildVertexArray(_lastFrame);

                // ****************
                // Handles all the OpenGL calls needed to display the point cloud
                // ****************

                // This clears the image, using what you set as GL.ClearColor earlier.
                // OpenGL provides several different types of data that can be rendered.
                // You can clear multiple buffers by using multiple bit flags.
                // However, we only modify the color, so ColorBufferBit is all we need to clear.
                GL.Clear(ClearBufferMask.ColorBufferBit);

                // Arguments:
                //   Which buffer the data should be sent to.
                //   How much data is being sent, in bytes. You can generally set this to the length of your array, multiplied by sizeof(array type).
                //   The vertices themselves.
                //   How the buffer will be used, so that OpenGL can write the data to the proper memory space on the GPU.
                //   There are three different BufferUsageHints for drawing:
                //     StaticDraw: This buffer will rarely, if ever, update after being initially uploaded.
                //     DynamicDraw: This buffer will change frequently after being initially uploaded.
                //     StreamDraw: This buffer will change on every frame.
                //   Writing to the proper memory space is important! Generally, you'll only want StaticDraw,
                //   but be sure to use the right one for your use case.
                GL.BufferData(BufferTarget.ArrayBuffer, _vertices.Length * sizeof(float), _vertices, BufferUsageHint.StaticDraw);

                // Bind the shader
                _shader.Use();

                // Bind the VAO
                GL.BindVertexArray(_vertexArrayObject);

                // Arguments:
                //   Primitive type; What sort of geometric primitive the vertices represent.
                //     OpenGL used to support many different primitive types, but almost all of the ones still supported
                //     is some variant of a triangle. Since we just want a single triangle, we use Triangles.
                //   Starting index; this is just the start of the data you want to draw. 0 here.
                //   How many vertices you want to draw. 3 for a triangle.
                GL.DrawArrays(PrimitiveType.Points, 0, _vertices.Length);
            }
        }

        private void GrabFrameHandler(object sender, RoutedEventArgs e)
        {
            GrabFrame();
        }


        bool _initialized = false;
        private void Init()
        {
            if (_initialized)
            {
                return;
            }

            GL.ClearColor(0.2f, 0.3f, 0.3f, 1.0f);
            _vertexBufferObject = GL.GenBuffer();
            GL.BindBuffer(BufferTarget.ArrayBuffer, _vertexBufferObject);
            _vertexArrayObject = GL.GenVertexArray();
            GL.BindVertexArray(_vertexArrayObject);

            _shader = new Shader("Shaders/shader.vert", "Shaders/shader.frag");
            _shader.Use();
            _camera = new Camera(Vector3.UnitZ * 3, 
                Convert.ToSingle(
                OpenTKControlWidth / OpenTKControlHeight));

            var vertexLocation = _shader.GetAttribLocation("aPosition");
            GL.EnableVertexAttribArray(vertexLocation);
            GL.VertexAttribPointer(vertexLocation, 3, VertexAttribPointerType.Float, false, 5 * sizeof(float), 0);

            var texCoordLocation = _shader.GetAttribLocation("aTexCoord");
            GL.EnableVertexAttribArray(texCoordLocation);
            GL.VertexAttribPointer(texCoordLocation, 2, VertexAttribPointerType.Float, false, 5 * sizeof(float), 3 * sizeof(float));

            _initialized = true;
        }

        private void OpenTkControl_Unloaded(object sender, RoutedEventArgs e)
        {

            // Unbind all the resources by binding the targets to 0/null.
            GL.BindBuffer(BufferTarget.ArrayBuffer, 0);
            GL.BindVertexArray(0);
            GL.UseProgram(0);

            // Delete all the resources.
            GL.DeleteBuffer(_vertexBufferObject);
            GL.DeleteVertexArray(_vertexArrayObject);

            GL.DeleteProgram(_shader.Handle);
        }

        private void MainWindow_KeyDown(object sender, System.Windows.Input.KeyEventArgs e)
        {
            const float cameraSpeed = 1.5f;

            if (e.Key == Key.W)
            {
                _camera.Position += _camera.Front * cameraSpeed; // Forward
            }

            if (e.Key == Key.S)
            {
                _camera.Position -= _camera.Front * cameraSpeed; // Backwards
            }
            if (e.Key == Key.A)
            {
                _camera.Position -= _camera.Right * cameraSpeed; // Left
            }
            if (e.Key == Key.D)
            {
                _camera.Position += _camera.Right * cameraSpeed; // Right
            }
            if (e.Key == Key.Space)
            {
                _camera.Position += _camera.Up * cameraSpeed; // Up
            }
            if (e.Key == Key.LeftShift)
            {
                _camera.Position -= _camera.Up * cameraSpeed; // Down
            }
        }

        private void StartStopCapture(object sender, RoutedEventArgs e)
        {
            if (_isCapturing)
            {
                _isCapturing = false;
                StartStop.Content = "Start";
            }
            else
            {
                _isCapturing = true;
                StartStop.Content = "Stop";
            }
        }

        System.Windows.Point _mouseStart;
        bool isMouseDown = false;

        private Point3DXY GetNearestPoint(double posX, double poxY)
        {
            if (_xyIndex == null)
            {
                return null;
            }

            double bestDis = double.MaxValue;
            var bestXS = -1f;
            var bestYS = -1f;
            Point3DXY result = null;
            foreach (var xd in _xyIndex)
            {
                foreach (var yd in xd.Value)
                {
                    foreach (var p in yd.Value)
                    {
                        var dx = p.SX - posX;
                        var dy = p.SY - poxY;
                        var dis2 = dx * dx + dy * dy;
                        if (dis2 < bestDis)
                        {
                            result = p;
                            bestXS = p.SX;
                            bestYS = p.SY;
                            bestDis = dis2;
                        }
                    }
                }
            }
            return result;
        }

        private double GetDistance(Point3D p1, Point3D p2)
        {
            var dx = p1.X - p2.X;
            var dy = p1.Y - p2.Y;
            var dz = p1.Z - p2.Z;
            return Math.Sqrt(dx * dx + dy * dy + dz * dz);
        }

        private void MainWindowMouseMove(object sender, MouseEventArgs e)
        {
            const float sensitivity = 0.2f;
            var currentPos = e.GetPosition(this);
            if (isMouseDown) // This bool variable is initially set to true.
            {
                var deltaX = currentPos.X - _mouseStart.X;
                var deltaY = currentPos.Y - _mouseStart.Y;
                _mouseStart = new System.Windows.Point(currentPos.X, currentPos.Y);
                bool isModelRotation = this.rotateModel.IsChecked ?? false;
                if (isModelRotation)
                {
                    _modelRotationX += Convert.ToSingle(deltaY)/100f;
                    _modelRotationY += Convert.ToSingle(deltaX)/100f;
                }
                else
                {
                    // Apply the camera pitch and yaw (we clamp the pitch in the camera class)
                    _camera.Yaw += Convert.ToSingle(deltaX * sensitivity);
                    _camera.Pitch -= Convert.ToSingle(deltaY * sensitivity); // Reversed since y-coordinates range from bottom to top
                }
            }

            if (!_isCapturing && _lastFrame != null)
            {
                _vertices = BuildVertexArray(_lastFrame);
                RefreshView();
            }

            if (_measure3DStart == null)
            {
                return;
            }

            var mousePosition = e.GetPosition(OpenTkControl);
            var nearest = GetNearestPoint(mousePosition.X, mousePosition.Y);
            if (nearest == null)
            {
                return;
            }

            var d = GetDistance(_measure3DStart, nearest);

            var measureDis = drawDistance.IsChecked ?? false;
            if (measureDis)
            {
                var canvasMousePos = e.GetPosition(drawCanvas);
                if (canvasMousePos.X >= 0
                    && canvasMousePos.Y >= 0
                    && canvasMousePos.X < drawCanvas.ActualWidth
                    && canvasMousePos.Y < drawCanvas.ActualHeight)
                {
                    distanceLine.X2 = canvasMousePos.X;
                    distanceLine.Y2 = canvasMousePos.Y;
                    measureText.Text = $"{d*100d:0.#}cm";

                    measureText.Foreground = new SolidColorBrush(Colors.Red);
                    Canvas.SetLeft(measureText, canvasMousePos.X - 10);
                    Canvas.SetTop(measureText, canvasMousePos.Y - 16);
                }
            }


            //_clipx = Convert.ToSingle(2.0 * mousePosition.X / OpenTkControl.ActualWidth - 1);
            //_clipy = Convert.ToSingle(2.0 * mousePosition.Y / OpenTkControl.ActualHeight - 1);

            //if (e.LeftButton == MouseButtonState.Pressed)
            //{
            //System.Windows.Point worldCoords = convertScreenToWorldCoords(Convert.ToInt32(currentPos.X), Convert.ToInt32(currentPos.Y));
            MousePosX.Text = $"X: {mousePosition.X},{nearest.SX}";
            MousePosY.Text = $"Y: {mousePosition.Y},{nearest.SY}";
            //}


        }

        private void MainWindowMouseUp(object sender, MouseButtonEventArgs e)
        {
            isMouseDown = false;
        }

        private void MainWindowMouseWheel(object sender, MouseWheelEventArgs e)
        {
            _camera.Fov -= e.Delta / 40;
        }

        private void GrabImageClicked(object sender, RoutedEventArgs e)
        {
            _grabImage = true;
            _grabDepthImage = true;
        }

        private void CutoutliersChanged(object sender, RoutedEventArgs e)
        {
            _cutOutliers = cutOutliers.IsChecked ?? false;
        }


        // functions:
        public System.Windows.Point convertScreenToWorldCoords(int x, int y)
        {
            int[] viewport = new int[4];
            //Matrix4 modelViewMatrix, projectionMatrix;
            //GL.GetFloat(GetPName.Modelview0MatrixExt, out modelViewMatrix); // GetPName.ModelviewMatrix
            //GL.GetFloat(GetPName.TransposeProjectionMatrix, out projectionMatrix); //GetPName.ProjectionMatrix

            Matrix4 modelViewMatrix = _camera.GetViewMatrix();
            Matrix4 projectionMatrix = _camera.GetProjectionMatrix();

            GL.GetInteger(GetPName.Viewport, viewport);
            Vector2 mouse;
            mouse.X = x;
            mouse.Y = y;

            Vector4 vector = UnProject(ref projectionMatrix, modelViewMatrix, new System.Windows.Size(viewport[2], viewport[3]), mouse);
            System.Windows.Point coords = new System.Windows.Point(vector.X, vector.Y);
            return coords;
        }

        public static Vector4 UnProject(ref Matrix4 projection, Matrix4 view, System.Windows.Size viewport, Vector2 mouse)
        {
            Vector4 vec;

            vec.X = 2.0f * mouse.X / (float)viewport.Width - 1;
            vec.Y = 2.0f * mouse.Y / (float)viewport.Height - 1;
            vec.Z = 0;
            vec.W = 1.0f;
            Matrix4 viewInv = Matrix4.Invert(view);
            Matrix4 projInv = Matrix4.Invert(projection);

            Vector4.Transform(ref vec, ref projInv, out vec);
            Vector4.Transform(ref vec, ref viewInv, out vec);

            if (vec.W > float.Epsilon || vec.W < float.Epsilon)
            {
                vec.X /= vec.W;
                vec.Y /= vec.W;
                vec.Z /= vec.W;
            }

            return vec;
        }

        private void GetScreenPos(
            float x, float y, float z,
            //float clipx, float clipy,
            float width, float height,
            out float xOut, out float yOut)
        {
            var model = Matrix4.Identity * Matrix4.CreateRotationX(0);
            var view = _camera.GetViewMatrix();
            var projection = _camera.GetProjectionMatrix();

            Vector4 vec;
            vec.X = x;
            vec.Y =  y;
            vec.Z = z;
            vec.W = 1.0f;

            Vector4.Transform(ref vec, ref model, out vec);
            Vector4.Transform(ref vec, ref view, out vec);
            Vector4.Transform(ref vec, ref projection, out vec);

            vec.X /= vec.W;
            vec.Y /= vec.W;

            xOut = Convert.ToSingle(((vec.X + 1.0) / 2.0) * width);
            yOut = Convert.ToSingle(((vec.Y + 1.0) / 2.0) * height);

            //var dxp = (vec.X - clipx) * width / 2f;
            //var dyp = (-vec.Y - clipy) * height / 2f;
            //return Math.Sqrt(dxp * dxp + dyp * dyp);
        }

        private void MainWindowMouseDown(object sender, MouseButtonEventArgs e)
        {
            var mousePosition = e.GetPosition(this);
            if (!isMouseDown)
            {
                _mouseStart = mousePosition;
                isMouseDown = true;
            }
        }

        private void OpenTkControl_SizeChanged(object sender, SizeChangedEventArgs e)
        {
            OpenTKControlWidth = OpenTkControl.ActualWidth;
            OpenTKControlHeight = OpenTkControl.ActualHeight;
        }

        private void GLWpfControlMouseLeftButtonDown(object sender, MouseButtonEventArgs e)
        {
            var measureDis = drawDistance.IsChecked ?? false;
            if (!measureDis)
            {
                return;
            }

            var mousePos = e.GetPosition(drawCanvas);
            _measure3DStart = GetNearestPoint(mousePos.X, mousePos.Y);
            if (_measure3DStart == null)
            {
                distanceLine.Visibility = Visibility.Collapsed;
                return;
            }

            distanceLine.Visibility = Visibility.Visible;
            //_measureStart = new System.Windows.Point(_measure3DStart.SX, _measure3DStart.SY);
            //distanceLine.X1 = _measureStart.X;
            //distanceLine.Y1 = _measureStart.Y;
            //distanceLine.X2 = _measureStart.X;
            //distanceLine.Y2 = _measureStart.Y;
        }

        private void UpdateMeasureStart()
        {

        }

        private void drawDistanceUnchecked(object sender, RoutedEventArgs e)
        {
            distanceLine.Visibility = Visibility.Collapsed;
        }

        private void UpdateRotCenterCheckedUnChecked(object sender, RoutedEventArgs e)
        {
            _updateRotCenter = updateRotCenter.IsChecked ?? false;
        }

        private void MainWindowClosed(object sender, EventArgs e)
        {
            if (_timer != null)
            {
                _timer.Stop();
                _timer.Dispose();
                _timer = null;
            }
        }
    }
}
