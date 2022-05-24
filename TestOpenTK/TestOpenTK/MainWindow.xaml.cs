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

namespace TestOpenTK
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        Timer _timer;
        bool _grabImage = false;
        int _nextImageIndex = 1;
        bool _cutOutliers = true;

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
        private float[] _vertices ={};

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

            Dispatcher.Invoke(new Action(() => {
                GrabFrame();
            }));

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
            Init();
            GL.Clear(ClearBufferMask.ColorBufferBit);
            GL.BufferData(BufferTarget.ArrayBuffer, _vertices.Length * sizeof(float), _vertices, BufferUsageHint.StaticDraw);

            if (_bmp != null)
            {
                _texture = Texture.LoadFromBitmap(_bmp);
                _texture.Use(TextureUnit.Texture0);
            }

            // Bind the shader
            _shader.Use();

            // Bind the VAO
            GL.BindVertexArray(_vertexArrayObject);

            var model = Matrix4.Identity * Matrix4.CreateRotationX(0);
            _shader.SetMatrix4("model", model);
            _shader.SetMatrix4("view", _camera.GetViewMatrix());
            _shader.SetMatrix4("projection", _camera.GetProjectionMatrix());

            GL.DrawArrays(PrimitiveType.Points, 0, _vertices.Length);
        }

        private float[] BuildVertexArray(
            float[] pointsX,
            float[] pointsY,
            float[] pointsZ,
            float[] textureU,
            float[] textureV)
        {
            List<float> res = new List<float>();
            for (int i = 0; i < pointsZ.Length; ++i)
            {
                if (pointsZ[i] != 0 && (!_cutOutliers ||
                    (textureU[i] >0 && textureU[i] <=1
                    && textureV[i] > 0 && textureV[i] <= 1)))
                {
                    res.Add(pointsX[i]);
                    res.Add(-pointsY[i]);
                    res.Add(pointsZ[i]);
                    res.Add(textureU[i]);
                    res.Add(textureV[i]);
                }
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

        private static void MySaveBMP(IntPtr intptr, int size, int width, int height)
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
            b.Save(@"c:\Misc\myPic.bmp", ImageFormat.Bmp);
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

            if (_bmp != null)
            {
                _bmp.Dispose();
                _bmp = null;
            }

            _bmp = BMPFromIntPtr(texturePixelsOut, textureWidthOut * textureHeightOut * 3, textureWidthOut, textureHeightOut);

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

            _vertices = BuildVertexArray(pointsX, pointsY, pointsZ, textureU, textureV);

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
            _camera = new Camera(Vector3.UnitZ * 3, Convert.ToSingle(OpenTkControl.ActualWidth / OpenTkControl.ActualHeight));

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

        private void MainWindowMouseMove(object sender, MouseEventArgs e)
        {
            const float sensitivity = 0.2f;
            var currentPos = e.GetPosition(this);
            if (isMouseDown) // This bool variable is initially set to true.
            {
                var deltaX = currentPos.X - _mouseStart.X;
                var deltaY = currentPos.Y - _mouseStart.Y;
                _mouseStart = new System.Windows.Point(currentPos.X, currentPos.Y);

                // Apply the camera pitch and yaw (we clamp the pitch in the camera class)
                _camera.Yaw += Convert.ToSingle(deltaX * sensitivity);
                _camera.Pitch -= Convert.ToSingle(deltaY * sensitivity); // Reversed since y-coordinates range from bottom to top
            }
        }

        private void MainWindowDown(object sender, MouseButtonEventArgs e)
        {
            if (isMouseDown)
            {
                return;
            }

            _mouseStart = e.GetPosition(this);
            isMouseDown = true;
        }

        private void MainWindowMouseUp(object sender, MouseButtonEventArgs e)
        {
            isMouseDown = false;
        }

        private void MainWindowMouseWheel(object sender, MouseWheelEventArgs e)
        {
            _camera.Fov -= e.Delta/40;
        }

        private void GrabImageClicked(object sender, RoutedEventArgs e)
        {
            _grabImage = true;
        }

        private void CutoutliersChanged(object sender, RoutedEventArgs e)
        {
            _cutOutliers = cutOutliers.IsChecked ?? false;
        }
    }
}
