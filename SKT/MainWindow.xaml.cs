using ScottPlot;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Media;
using System.Windows.Forms;
using MessageBox = System.Windows.MessageBox;

namespace SKT
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();

            LabelCN.DataContext = mesh_forward;
            LabelMinX.DataContext = mesh_forward;
            LabelMaxX.DataContext = mesh_forward;
            LabelMinZ.DataContext = mesh_forward;
            LabelMaxZ.DataContext = mesh_forward;

            LabelCE.DataContext = mesh_forward;
            LabelMes.DataContext = mesh_forward;

            LabelMinP.DataContext = mesh_forward;
            LabelMaxP.DataContext = mesh_forward;

            MeshChart.LeftClicked += MeshChart_LeftClicked;

            TextBoxStartGamma.DataContext = inverse_input;
            TextBoxAllowMaxDifference.DataContext = inverse_input;
            TextBoxGammaMul.DataContext = inverse_input;
            TextBoxMaxFunctionalInPercent.DataContext = inverse_input;
            TextBoxMaxIterations.DataContext = inverse_input;

            ListOfInverseIteration.DataContext = inverse_results;
        }

        RectangleMesh mesh_forward = new();
        RectangleMesh mesh_inverse = new();

        List<Profile_receivers> profiles_forward = new();
        List<Profile_receivers> profiles_source = new();
        List<Profile_receivers> profiles_syntetic = new();
        List<Profile_receivers> profiles_difference = new();

        InverseInput inverse_input = new();
        InverseResult inverse_results = new();

        private void Read_mesh_from_file(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog = new Microsoft.Win32.OpenFileDialog();
            if (openFileDialog.ShowDialog() == true)
            {
                mesh_forward.Clear();
                var s = File.ReadAllLines(openFileDialog.FileName);

                var coords = s[0].Split('\t');
                int count_x = int.Parse(coords[0]);
                int count_z = int.Parse(coords[1]);

                coords = s[1].Split('\t');
                double min_x = double.Parse(coords[0]);
                double max_x = double.Parse(coords[1]);

                coords = s[2].Split('\t');
                double min_z = double.Parse(coords[0]);
                double max_z = double.Parse(coords[1]);

                mesh_forward.BuildMesh(count_x, count_z, min_x, max_x, min_z, max_z);
            }
        }

        private void Generate_mesh(object sender, RoutedEventArgs e)
        {
            GenerateMeshWindow window = new GenerateMeshWindow();
            if (window.ShowDialog() == true)
            {
                mesh_forward.Clear();
                mesh_forward.BuildMesh(window.NX, window.NZ, window.MinX, window.MaxX, window.MinZ, window.MaxZ);

                NodesGrid.ItemsSource = mesh_forward.GetPoints();
                ElementsGrid.ItemsSource = mesh_forward.GetElements();
                PGrid.ItemsSource = mesh_forward.P;
                
                foreach (DoubleValue p in mesh_forward.P)
                {
                    p.PropertyChanged += P_PropertyChanged;
                }

                MeshChart.Plot.Clear();
                DrawMediumProperties(MeshChart, mesh_forward);
                DrawMesh(MeshChart, mesh_forward);
                MeshChart.Plot.AxisAuto();
                MeshChart.Plot.XAxis.Label("X");
                MeshChart.Plot.YAxis.Label("Z");
                MeshChart.Refresh();

                UpdateForwardModel(save_axis_limits: false);
            }
        }

        private void MeshChart_LeftClicked(object sender, RoutedEventArgs e)
        {
            if ((mesh_forward.X.Count > 1) && (mesh_forward.Z.Count > 1))
            {
                (double x, double z) = MeshChart.GetMouseCoordinates();

                if ((x >= mesh_forward.MinX) && (x <= mesh_forward.MaxX) && (z >= mesh_forward.MinZ) && (z <= mesh_forward.MaxZ))
                {
                    double h_x = (mesh_forward.MaxX - mesh_forward.MinX) / (mesh_forward.X.Count - 1);
                    double h_z = (mesh_forward.MaxZ - mesh_forward.MinZ) / (mesh_forward.Z.Count - 1);

                    int index_i = (int)((x - mesh_forward.MinX) / h_x);
                    int index_j = (int)((z - mesh_forward.MinZ) / h_z);

                    int index_p = index_i * (mesh_forward.Z.Count - 1) + index_j;

                    var cell = new DataGridCellInfo(PGrid.Items[index_p], PGrid.Columns[1]);

                    PGrid.SelectedIndex = index_p;
                    PGrid.CurrentCell = cell;

                    PGrid.Focus();
                    PGrid.ScrollIntoView(PGrid.SelectedItem);
                    PGrid.BeginEdit(e);
                }
            }
        }

        private void P_PropertyChanged(object? sender, System.ComponentModel.PropertyChangedEventArgs e)
        {
            var axis = MeshChart.Plot.GetAxisLimits();

            MeshChart.Plot.Clear();

            DrawMediumProperties(MeshChart, mesh_forward);
            DrawMesh(MeshChart, mesh_forward);

            MeshChart.Plot.SetAxisLimits(axis);
            MeshChart.Plot.XAxis.Label("X");
            MeshChart.Plot.YAxis.Label("Z");
            MeshChart.Refresh();

            UpdateForwardModel(save_axis_limits: true);
        }

        private void Save_mesh(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.SaveFileDialog saveFileDialog = new Microsoft.Win32.SaveFileDialog();
            saveFileDialog.Filter = "Text files (*.txt)|*.txt";
            if (saveFileDialog.ShowDialog() == true)
            {
                var file = new StreamWriter(saveFileDialog.FileName);
                file.WriteLine($"{mesh_forward.X.Count}\t{mesh_forward.Z.Count}");
                file.WriteLine($"{mesh_forward.MinX}\t{mesh_forward.MaxX}");
                file.WriteLine($"{mesh_forward.MinZ}\t{mesh_forward.MaxZ}");
                file.Close();
            }
        }


        private void Read_receivers_from_file(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.OpenFileDialog openFileDialog = new Microsoft.Win32.OpenFileDialog();
            if (openFileDialog.ShowDialog() == true)
            {
                mesh_forward.Clear();
                var s = File.ReadAllLines(openFileDialog.FileName);

                int count_nodes = int.Parse(s[0]);

                for (int i = 1; i < count_nodes + 1; i++)
                {
                    var coords = s[i].Split('\t');
                    //mesh.Nodes.Add(new Point(double.Parse(coords[0]), double.Parse(coords[1])));
                }

                int count_elements = int.Parse(s[count_nodes]);

                for (int i = count_nodes; i < count_elements + count_nodes + 1; i++)
                {
                    var str = s[i].Split('\t');
                    //mesh.Elements.Add(new Element(
                    //    int.Parse(str[0]),
                    //    int.Parse(str[1]),
                    //    int.Parse(str[2]),
                    //    int.Parse(str[3])));
                }
            }
        }

        private void Generate_receivers(object sender, RoutedEventArgs e)
        {
            GenerateProfilesWindow window = new GenerateProfilesWindow();
            if (window.ShowDialog() == true)
            {
                double x_0 = window.X0;
                double z_0 = window.Z0;
                double h_x = window.HX;
                double h_z = window.HZ;
                double offset_x = window.OffsetX;
                double offset_z = window.OffsetZ;
                int count_point = window.CountPoint;

                foreach (var profile in profiles_forward)
                {
                    profile.Values.Clear();
                    profile.Coords.Clear();
                }

                NodesOfProfilesGrid.ItemsSource = null;
                ProfilesGrid.ItemsSource = null;

                profiles_forward.Clear();
                ObservingChart.Plot.Clear();

                int num = 1;

                for (int i_nx = 0; i_nx < window.NX; i_nx++)
                {
                    for (int j_nz = 0; j_nz < window.NZ; j_nz++, num++)
                    {
                        Profile_receivers profile = new()
                        {
                            X0 = offset_x * i_nx + x_0,
                            Z0 = offset_z * j_nz + z_0,
                            X1 = offset_x * i_nx + x_0 + h_x,
                            Z1 = offset_z * j_nz + z_0 + h_z,
                            N = num,

                            Count = count_point
                        };

                        profile.calculate_coords();

                        double[] xs = new double[profile.Coords.Count];
                        double[] zs = new double[profile.Coords.Count];

                        for (int i = 0; i < profile.Coords.Count; i++)
                        {
                            xs[i] = profile.Coords[i].X;
                            zs[i] = profile.Coords[i].Z;
                        }

                        var color = ObservingChart.Plot.GetNextColor();
                        profile.Brush = new SolidColorBrush(System.Windows.Media.Color.FromArgb(color.A, color.R, color.G, color.B));
                        var scatter = ObservingChart.Plot.AddScatter(xs, zs, color);

                        profiles_forward.Add(profile);
                    }
                }

                ProfilesGrid.ItemsSource = profiles_forward;
                ObservingChart.Plot.XAxis.Label("X");
                ObservingChart.Plot.YAxis.Label("Z");
                ObservingChart.Refresh();

                UpdateForwardModel(save_axis_limits: false);
            }

        }

        private void Save_receivers(object sender, RoutedEventArgs e)
        {
            Microsoft.Win32.SaveFileDialog saveFileDialog = new Microsoft.Win32.SaveFileDialog();
            saveFileDialog.Filter = "Text files (*.txt)|*.txt";
            if (saveFileDialog.ShowDialog() == true)
            {
                var file = new StreamWriter(saveFileDialog.FileName);
                //file.WriteLine($"{mesh.Nodes.Count}");

                //foreach (Point point in mesh.Nodes)
                //{
                //    file.WriteLine($"{point.X}\t{point.Z}");
                //}

                //file.WriteLine($"{mesh.Elements.Count}");

                //foreach (Element num_nodes in mesh.Elements)
                //{
                //    file.WriteLine($"{num_nodes.Node1}\t{num_nodes.Node2}\t{num_nodes.Node3}\t{num_nodes.Node4}");
                //}

                file.Close();
            }
        }

        private void ProfilesGrid_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            int selected_index = ProfilesGrid.SelectedIndex;
            if (selected_index >= 0 && selected_index < profiles_forward.Count)
            {
                NodesOfProfilesGrid.ItemsSource = profiles_forward[selected_index].Coords;
            }
        }

        private void CheckBoxShowObserving_Checked(object sender, RoutedEventArgs e)
        {
            UpdateForwardModel(save_axis_limits: true);
        }

        private void CheckBoxShowObserving_Unchecked(object sender, RoutedEventArgs e)
        {
            UpdateForwardModel(save_axis_limits: true);
        }

        private void CheckBoxShowModel_Checked(object sender, RoutedEventArgs e)
        {
            UpdateForwardModel(save_axis_limits: true);
        }

        private void CheckBoxShowModel_Unchecked(object sender, RoutedEventArgs e)
        {
            UpdateForwardModel(save_axis_limits: true);
        }

        private void Forward_Calculation(object sender, RoutedEventArgs e)
        {
            if (ComboBoxCalculationMode.SelectedItem != null)
            {
                Solver.SolveForward(mesh_forward, profiles_forward, (CalculationMode)ComboBoxCalculationMode.SelectedItem);

                UpdateForwardChart(save_axis_limits: false);
            }
            else
            {
                MessageBox.Show("Select calculation mode!");
            }
        }

        private void UpdateObservingChart(bool save_axis_limits)
        {
            var axis = ObservingChart.Plot.GetAxisLimits();

            ObservingChart.Plot.Clear();

            DrawProfiles(ObservingChart, profiles_forward);

            if (save_axis_limits)
            {
                ObservingChart.Plot.SetAxisLimits(axis);
            }
            else
            {
                ObservingChart.Plot.AxisAuto();
            }

            ObservingChart.Plot.XAxis.Label("X");
            ObservingChart.Plot.YAxis.Label("Z");

            ObservingChart.Refresh();
        }

        private void UpdateForwardModel(bool save_axis_limits)
        {
            var axis = ForwardModel.Plot.GetAxisLimits();
            ForwardModel.Plot.Clear();

            if (CheckBoxShowModel.IsChecked == true)
            {
                DrawMediumProperties(ForwardModel, mesh_forward);
                DrawMesh(ForwardModel, mesh_forward);
            }

            if (CheckBoxShowObserving.IsChecked == true)
            {
                DrawProfiles(ForwardModel, profiles_forward);
            }

            if (save_axis_limits)
            {
                ForwardModel.Plot.SetAxisLimits(axis);
            }
            else
            {
                ForwardModel.Plot.AxisAuto();
            }

            ForwardModel.Plot.XAxis.Label("X");
            ForwardModel.Plot.YAxis.Label("Z");

            ForwardModel.Refresh();
        }

        private void UpdateInverseModelChart(bool save_axis_limits)
        {
            var axis = InverseModelChart.Plot.GetAxisLimits();

            InverseModelChart.Plot.Clear();

            DrawMediumProperties(InverseModelChart, mesh_inverse);
            DrawMesh(InverseModelChart, mesh_inverse);
            DrawProfiles(InverseModelChart, profiles_source);

            if (save_axis_limits)
            {
                InverseModelChart.Plot.SetAxisLimits(axis);
            }
            else
            {
                InverseModelChart.Plot.AxisAuto();
            }

            InverseModelChart.Plot.XAxis.Label("X");
            InverseModelChart.Plot.YAxis.Label("Z");

            InverseModelChart.Refresh();
        }

        private void UpdateForwardChart(bool save_axis_limits)
        {
            var axis = ForwardChart.Plot.GetAxisLimits();

            ForwardChart.Plot.Clear();

            DrawValuesInProfiles(ForwardChart, profiles_forward);

            if (save_axis_limits)
            {
                ForwardChart.Plot.SetAxisLimits(axis);
            }
            else
            {
                ForwardChart.Plot.AxisAuto();
            }
            
            ForwardChart.Plot.XAxis.Label("X, m");
            ForwardChart.Plot.YAxis.Label("Bx, T");

            ForwardChart.Refresh();
        }

        private void ButtonChangeColor_Click(object sender, RoutedEventArgs e)
        {
           ColorDialog MyDialog = new ColorDialog();
           if (MyDialog.ShowDialog() == System.Windows.Forms.DialogResult.OK)
           {
                var color = MyDialog.Color;

                var b = sender as System.Windows.Controls.Button;
                if (b != null)
                {
                    b.Background = new SolidColorBrush(System.Windows.Media.Color.FromArgb(color.A, color.R, color.G, color.B));
                    UpdateObservingChart(save_axis_limits: true);
                    UpdateForwardModel(save_axis_limits: true);
                    UpdateForwardChart(save_axis_limits: true);
                }
            }
        }

        private void LoadFromObserving(object sender, RoutedEventArgs e)
        {
            if (profiles_forward.Count == 0)
            {
                MessageBox.Show("Forward result is not loaded or calculated!");
                return;
            }

            profiles_source.Clear();
            profiles_syntetic.Clear();
            profiles_difference.Clear();

            for (int i = 0; i < profiles_forward.Count; i++)
            {
                Profile_receivers profile_1 = new();
                profile_1.DeepCopy(profiles_forward[i]);
                profiles_source.Add(profile_1);

                Profile_receivers profile_2 = new();
                profile_2.Copy(profiles_forward[i]);
                profiles_syntetic.Add(profile_2);

                Profile_receivers profile_3 = new();
                profile_3.Copy(profiles_forward[i]);
                profiles_difference.Add(profile_3);
            }

            UpdateInverseModelChart(save_axis_limits: true);

            ReadDataChart.Plot.Clear();
            DrawValuesInProfiles(ReadDataChart, profiles_source, "Input data");
            ReadDataChart.Plot.XAxis.Label("X");
            ReadDataChart.Plot.YAxis.Label("Bx");
            ReadDataChart.Refresh();

            InverseDataChart.Plot.Clear();
            InverseDataChart.Refresh();

            DeltaDataChart.Plot.Clear();
            DeltaDataChart.Refresh();
        }

        private void LoadFromModel(object sender, RoutedEventArgs e)
        {
            if (mesh_forward.CountElements <= 0)
            {
                MessageBox.Show("Mesh is not loaded or generated!");
                return;
            }

            mesh_inverse.Clear();
            mesh_inverse.BuildMesh(mesh_forward.X.Count, mesh_forward.Z.Count, mesh_forward.MinX, mesh_forward.MaxX, mesh_forward.MinZ, mesh_forward.MaxZ);

            UpdateInverseModelChart(save_axis_limits: false);
        }

        public void DrawMesh(WpfPlot plotter, RectangleMesh mesh)
        {
            List<List<(double x, double y)>> polys = new List<List<(double x, double y)>>();

            for (int i_x = 0; i_x < mesh.X.Count - 1; i_x++)
            {
                for (int j_z = 0; j_z < mesh.Z.Count - 1; j_z++)
                {
                    double[] xs = { mesh.X[i_x], mesh.X[i_x + 1], mesh.X[i_x + 1], mesh.X[i_x] };
                    double[] ys = { mesh.Z[j_z + 1], mesh.Z[j_z + 1], mesh.Z[j_z], mesh.Z[j_z] };

                    List<(double x, double y)> thisPolygon = xs.Zip(ys, (xp, yp) => (xp, yp)).ToList();
                    polys.Add(thisPolygon);
                }
            }

            plotter.Plot.AddPolygons(polys, plotter.Plot.GetNextColor(0.0), lineWidth: 1, System.Drawing.Color.Gray);

            for (int i_x = 0; i_x < mesh.X.Count - 1; i_x++)
            {
                for (int j_z = 0; j_z < mesh.Z.Count - 1; j_z++)
                {
                    double x = (mesh.X[i_x]  + mesh.X[i_x + 1]) / 2.0;
                    double y = (mesh.Z[j_z] + mesh.Z[j_z + 1]) / 2.0;

                    double diff;
                    if (mesh.MaxP != 0.0)
                    {
                        diff = mesh.P[i_x * (mesh.Z.Count - 1) + j_z].Value / mesh.MaxP;
                        
                    }
                    else
                    {
                        diff = 1;
                    }

                    if(diff > mesh_forward.MaxP)
                    {
                        diff = 1;
                    }

                    if (diff < 0)
                    {
                        diff = 0;
                    }

                    var labelFont = new ScottPlot.Drawing.Font
                    {
                        Bold = true,
                        Color = System.Drawing.Color.FromArgb((int)(255 * diff), (int)(255 * diff), (int)(255 * diff)),
                        Alignment = Alignment.MiddleCenter
                    };

                    plotter.Plot.AddText(mesh.P[i_x * (mesh.Z.Count - 1) + j_z].Value.ToString("F2"), x, y, labelFont);
                }
            }
        }

        public void DrawMediumProperties(WpfPlot plotter, RectangleMesh mesh)
        {
            if ((mesh.CountNodes != 0) && (mesh.MaxP > 0))
            {
                double?[,] medium_properties_data = new double?[mesh.Z.Count - 1, mesh.X.Count - 1];

                for (int i = 0; i < mesh.X.Count - 1; i++)
                {
                    for (int j = 0; j < mesh.Z.Count - 1; j++)
                    {
                        medium_properties_data[j, i] = mesh.P[i * (mesh.Z.Count - 1) + j].Value;
                    }
                }

                var hm = plotter.Plot.AddHeatmap(medium_properties_data, ScottPlot.Drawing.Colormap.GrayscaleR);
                hm.Update(medium_properties_data, min: 0, max: mesh.MaxP);
                var cb = plotter.Plot.AddColorbar(hm);

                cb.MinValue = 0;
                cb.MaxValue = mesh_forward.MaxP;

                hm.FlipVertically = true;
                hm.OffsetX = mesh.MinX;
                hm.OffsetY = mesh.MinZ;

                hm.CellWidth = (mesh.MaxX - mesh.MinX) / (mesh.X.Count - 1.0);
                hm.CellHeight = (mesh.MaxZ - mesh.MinZ) / (mesh.Z.Count - 1.0);
            }
        }

        public void DrawProfiles(WpfPlot plotter, List<Profile_receivers> profiles)
        {
            foreach (Profile_receivers profile in profiles)
            {
                double[] xs = new double[profile.Coords.Count];
                double[] zs = new double[profile.Coords.Count];

                for (int i = 0; i < profile.Coords.Count; i++)
                {
                    xs[i] = profile.Coords[i].X;
                    zs[i] = profile.Coords[i].Z;
                }

                var color = profile.Brush.Color;
                plotter.Plot.AddScatter(xs, zs, System.Drawing.Color.FromArgb(color.A, color.R, color.G, color.B));
            }
        }

        private void DrawValuesInProfiles(WpfPlot plotter, List<Profile_receivers> profiles, string title="")
        {
            foreach (var profile in profiles)
            {
                double[] xs = new double[profile.Count];
                double[] zs = new double[profile.Count];

                for (int i = 0; i < profile.Count; i++)
                {
                    xs[i] = profile.Coords[i].X;
                    zs[i] = profile.Values[i];
                }
                var color = profile.Brush.Color;
                plotter.Plot.AddScatter(xs, zs, System.Drawing.Color.FromArgb(color.A, color.R, color.G, color.B));
            }
            plotter.Plot.Title(title);
        }

        private void Inverse_Calculation(object sender, RoutedEventArgs e)
        {
            if (mesh_inverse.CountElements < 1)
            {
                MessageBox.Show("Mesh is not loaded or generated");
                return;
            }

            if (profiles_source.Count == 0)
            {
                MessageBox.Show("Observing is not loaded or generated!");
                return;
            }

            if (ComboBoxInverseCalculationMode.SelectedItem == null)
            {
                MessageBox.Show("Select calculation mode!");
                return;
            }

            if (CheckBoxRegularization.IsChecked == true)
            {
                inverse_results.Clear();

                inverse_results = Solver.SolveInverse(
                    mesh_inverse,
                    profiles_source,
                    (CalculationMode)ComboBoxInverseCalculationMode.SelectedItem,
                    inverse_input.StartGamma,
                    inverse_input.AllowMaxDifference,
                    inverse_input.GammaMult,
                    inverse_input.MaxFunctionalInPercent,
                    inverse_input.MaxIterations,
                    inverse_input.SaveAllIterations
                );
            }
            else
            {
                inverse_results = Solver.SolveInverse(
                    mesh_inverse,
                    profiles_source,
                    (CalculationMode)ComboBoxInverseCalculationMode.SelectedItem,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1,
                    true
                );
            }

            ListOfInverseIteration.ItemsSource = inverse_results.Iters;

            UpdateInverseModelChart(save_axis_limits: true);
        }

        private void RadioButton_ShowAllIterations(object sender, RoutedEventArgs e)
        {
            inverse_input.SaveAllIterations = true;
        }

        private void RadioButton_ShowFirstAndLast(object sender, RoutedEventArgs e)
        {
            inverse_input.SaveAllIterations = false;
        }

        private void ListOfInverseIteration_SelectionChanged(object sender, SelectionChangedEventArgs e)
        {
            int selected_index = ListOfInverseIteration.SelectedIndex;
            if (selected_index >= 0 && selected_index < inverse_results.Iters.Count)
            {
                ListOfInverseP.ItemsSource = inverse_results.Iters[selected_index].ValuesP;

                var p = inverse_results.Iters[selected_index].ValuesP;
                for (int i = 0; i < mesh_inverse.CountElements; i++)
                {
                    mesh_inverse.P[i].Value = p[i];
                }
                //mesh_inverse.MinP = 0;
                //mesh_inverse.MaxP = mesh_forward.MaxP;

                var b = inverse_results.Iters[selected_index].ValuesB;
                for (int i = 0; i < profiles_syntetic.Count; i++)
                {
                    var profile = profiles_syntetic[i];
                    var profile_b = b[i];
                    for (int j = 0; j < profile.Count; j++)
                    {
                        profile.Values[j] = profile_b[j];
                    }
                }

                for (int i = 0; i < profiles_syntetic.Count; i++)
                {
                    var profile_syntetic = profiles_syntetic[i];
                    var profile_source = profiles_source[i];
                    var profile_difference = profiles_difference[i];
                    for (int j = 0; j < profile_syntetic.Count; j++)
                    {
                        profile_difference.Values[j] = profile_source.Values[j] - profile_syntetic.Values[j];
                    }
                }

                UpdateInverseModelChart(save_axis_limits: true);

                InverseDataChart.Plot.Clear();
                DrawValuesInProfiles(InverseDataChart, profiles_syntetic, "Syntetic data");
                InverseDataChart.Refresh();

                DeltaDataChart.Plot.Clear();
                DrawValuesInProfiles(DeltaDataChart, profiles_difference, "Difference between input and syntetic data");
                DeltaDataChart.Refresh();
            }
        }

        void DataGrid_LoadingRow(object sender, DataGridRowEventArgs e)
        {
            e.Row.Header = (e.Row.GetIndex() + 1).ToString();
        }

        private void Generate_inverse_mesh(object sender, RoutedEventArgs e)
        {
            GenerateMeshWindow window = new();
            if (window.ShowDialog() == true)
            {
                mesh_inverse.Clear();
                mesh_inverse.BuildMesh(window.NX, window.NZ, window.MinX, window.MaxX, window.MinZ, window.MaxZ);

                UpdateInverseModelChart(save_axis_limits: false);
            }
        }
    }
}
