using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.ComponentModel;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows.Media;
using System.Windows.Media.Media3D;

namespace SKT
{
    public class Profile_receivers
    {
        public Profile_receivers() 
        {
            Coords = new();
            Values = new();

            count_ = 0;
        }

        public int N { get; set; }

        public System.Windows.Media.SolidColorBrush Brush { get; set; }

        // список координат точек измерения
        public List<Point> Coords { get; }
        
        // значения поля в точках измерения
        public List<double> Values { get; }

        private int count_ = 0;

        private double x0_ = 0;
        private double z0_ = 0;
        private double x1_ = 0;
        private double z1_ = 0;

        public double X0 { 
            get
            {
                return x0_;
            }
            set
            {
                x0_ = value;
            }
        }

        public double X1
        {
            get
            {
                return x1_;
            }
            set
            {
                x1_ = value;
            }
        }

        public double Z0
        {
            get
            {
                return z0_;
            }
            set
            {
                z0_ = value;
            }
        }

        public double Z1
        {
            get
            {
                return z1_;
            }
            set
            {
                z1_ = value;
            }
        }

        public int Count 
        { 
            get
            {
                return count_;
            }
            set
            {
                if (value != count_ && value > 0)
                {
                    count_ = value;
                }
            }
        }

        public void Clear()
        {
            Coords.Clear();
            Values.Clear();

            count_ = 0;
            N = 0;

            x0_ = 0;
            z0_ = 0;
            x1_ = 0;
            z1_ = 0;
        }

        public void DeepCopy(Profile_receivers profile)
        {
            this.Clear();

            Count = profile.Count;
            X0 = profile.X0;
            Z0 = profile.Z0;
            Z1 = profile.Z1;
            X1 = profile.X1;
            N = profile.N;

            Brush = profile.Brush;

            calculate_coords();

            for (int i = 0; i < profile.Values.Count; i++)
            {
                Values[i] = profile.Values[i];
            }
        }

        public void Copy(Profile_receivers profile)
        {
            this.Clear();

            Count = profile.Count;
            X0 = profile.X0;
            Z0 = profile.Z0;
            Z1 = profile.Z1;
            X1 = profile.X1;
            N = profile.N;

            Brush = profile.Brush;

            calculate_coords();

            for (int i = 0; i < profile.Values.Count; i++)
            {
                Values[i] = 0.0;
            }
        }

        public void calculate_coords()
        {
            Coords.Clear(); 
            Values.Clear();

            if ((X0 == X1) && (Z0 == Z1))
            {
                count_ = 1;
                Coords.Add(new Point(X0, Z0, 1));
                Values.Add(0.0);
                return;
            }

            double stepX = (X1 - X0) / (count_ - 1.0);
            double stepZ = (Z1 - Z0) / (count_ - 1.0);

            for (int i = 0; i < count_; i++)
            {
                Coords.Add(new Point(stepX * i + X0, stepZ * i + Z0, i + 1));
                Values.Add(0.0);
            }
        }
    }
}
