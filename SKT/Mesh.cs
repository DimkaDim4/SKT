using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
using System.Collections.Specialized;
using System.ComponentModel;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Media;
using System.Windows.Media.Media3D;

namespace SKT
{
    public class DoubleValue : INotifyPropertyChanged
    {
        double value_;

        public double Value {
            get { return value_; }
            set { 
                value_ = value;
                OnPropertyChanged("Value");
            }
        }

        public int N { get; set; }

        public event PropertyChangedEventHandler? PropertyChanged;

        public void OnPropertyChanged([CallerMemberName] string prop = "")
        {
            if (PropertyChanged != null)
                PropertyChanged(this, new PropertyChangedEventArgs(prop));
        }
    }

    public class Point
    {
        public int N { get; set; }
        public double X { get; set; }
        public double Z { get; set; }

        public Point()
        {
            this.X = 0;
            this.Z = 0;
            this.N = 0;
        }

        public Point(double x, double z, int n)
        {
            this.X = x;
            this.Z = z;
            this.N = n;
        }
    }

    public class Element
    {
        public int N { get; set; }
        public int N1 { get; set; }
        public int N2 { get; set; }
        public int N3 { get; set; }
        public int N4 { get; set; }

        public Element()
        {
            N = 0;
            N1 = 0;
            N2 = 0;
            N3 = 0;
            N4 = 0;
        }

        public Element(int n, int n_1, int n_2, int n_3, int n_4)
        {
            N = n;
            N1 = n_1;
            N2 = n_2;
            N3 = n_3;
            N4 = n_4;
        }
    }

    public class RectangleMesh : INotifyPropertyChanged
    {
        public RectangleMesh() 
        {
            X = new();
            Z = new();
            P = new();
            P.CollectionChanged += P_CollectionChanged;

            CountNodes = 0;
            CountElements = 0;
        }

        // Список узлов
        public List<double> X { get; private set; }
        public List<double> Z { get; private set; }
        public ObservableCollection<DoubleValue> P { get; set; }

        private void P_CollectionChanged(object? sender, System.Collections.Specialized.NotifyCollectionChangedEventArgs e)
        {
            if (e.Action == NotifyCollectionChangedAction.Add)
            {
                if (e.NewItems?[0] is DoubleValue value_)
                {
                    value_.PropertyChanged += P_PropertyChanged;

                    MinP = MinP < value_.Value ? MinP : value_.Value;
                    MaxP = MaxP > value_.Value ? MaxP : value_.Value;

                    OnPropertyChanged("MinP");
                    OnPropertyChanged("MaxP");
                }
            }
        }

        private void P_PropertyChanged(object? sender, PropertyChangedEventArgs e)
        {
            double new_value = (sender as DoubleValue).Value;

            MinP = P[0].Value;
            MaxP = P[0].Value;

            foreach (var p in P) 
            {
                MinP = MinP < p.Value ? MinP : p.Value;
                MaxP = MaxP > p.Value ? MaxP : p.Value;
            }

            OnPropertyChanged("MinP");
            OnPropertyChanged("MaxP");
        }

        public double MinP { get; set; }
        public double MaxP { get; set; }

        public double MinX { get; private set; }
        public double MaxX { get; private set; }

        public double MinZ { get; private set; }
        public double MaxZ { get; private set; }

        public double Mes { get; private set; }

        public int CountNodes { get; private set; }

        public int CountElements { get; private set; }

        public event PropertyChangedEventHandler? PropertyChanged;

        public void OnPropertyChanged([CallerMemberName] string prop = "")
        {
            if (PropertyChanged != null)
                PropertyChanged(this, new PropertyChangedEventArgs(prop));
        }

        public void Clear()
        {
            X.Clear();
            Z.Clear();
            P.Clear();

            MinP = 0;
            MaxP = 0;

            MinX = 0;
            MinZ = 0;
            MaxX = 0;
            MaxZ = 0;

            CountElements = 0;
            CountNodes = 0;

            Mes = 0;
        }

        public void BuildMesh(int Nx, int Nz, double MinX, double MaxX, double MinZ, double MaxZ)
        {
            X.Clear();
            Z.Clear();
            P.Clear();

            double h_x = (MaxX - MinX) / (Nx - 1.0);
            double h_z = (MaxZ - MinZ) / (Nz - 1.0);

            for (int i = 0; i < Nx; i++)
            {
                X.Add(i * h_x + MinX);
            }

            for (int i = 0; i < Nz; i++)
            {
                Z.Add(i * h_z + MinZ);
            }

            for (int i = 0; i < (Nx - 1) * (Nz - 1); i++)
            {
                P.Add(new DoubleValue { Value = 0, N = i + 1 });
            }

            Mes = h_x * h_z;

            this.MinX = X[0];
            this.MinZ = Z[0];
            this.MaxX = X[^1];
            this.MaxZ = Z[^1];

            CountNodes = X.Count * Z.Count;
            CountElements = (X.Count - 1) * (Z.Count - 1);

            OnPropertyChanged("MinX");
            OnPropertyChanged("MaxX");
            OnPropertyChanged("MinZ");
            OnPropertyChanged("MaxZ");
            OnPropertyChanged("CountNodes");
            OnPropertyChanged("CountElements");
            OnPropertyChanged("Mes");
        }

        public List<Point> GetPoints()
        {
            List<Point> points = new();

            for (int i = 0; i < X.Count; i++)
            {
                for (int j = 0; j < Z.Count; j++)
                {
                    points.Add( new Point(X[i], Z[j], i * Z.Count + j + 1));
                }
            }

            return points;
        }

        public List<Element> GetElements()
        {
            List<Element> elements = new();

            int nz = Z.Count;

            for (int i = 0; i < X.Count - 1; i++)
            {
                for (int j = 0; j < Z.Count - 1; j++)
                {
                    elements.Add(new Element(
                        i * (Z.Count - 1) + j + 1,
                        i * nz + j + 2,
                        (i + 1) * nz + j + 2,
                        (i + 1) * nz + j + 1,
                        i * nz + j + 1));
                }
            }

            return elements;
        }
    }
}
