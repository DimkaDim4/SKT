using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;

namespace SKT
{
    /// <summary>
    /// Логика взаимодействия для GenerateProfilesWindow.xaml
    /// </summary>
    public partial class GenerateProfilesWindow : Window
    {
        public GenerateProfilesWindow()
        {
            InitializeComponent();
            DataContext = this;
        }

        int nx = 1, nz = 1, count_point = 500;
        double x_0 = -1000, z_0 = 0, h_x = 2000, h_z = 0, offset_x = 0, offset_z = 0;

        public List<Profile_receivers> profiles = new();

        public int NX { get { return nx; } set { if (value > 0) nx = value; } }
        public int NZ { get { return nz; } set { if (value > 0) nz = value; } }

        public int CountPoint { get { return count_point; } set { if (value > 0) count_point = value; } }

        public double X0 { get { return x_0; } set { x_0 = value; } }

        public double Z0 { get { return z_0; } set { z_0 = value; } }

        public double HX { get { return h_x; } set { if (value >= 0) h_x = value; } }

        public double HZ { get { return h_z; } set { if (value >= 0) h_z = value; } }

        public double OffsetX { get { return offset_x; } set { offset_x = value; } }

        public double OffsetZ { get { return offset_z; } set { offset_z = value; } }



        private void DOUBLE_TextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            var textBox = sender as TextBox;
            double var;

            if (double.TryParse(textBox.Text, NumberStyles.Number, CultureInfo.InvariantCulture, out var))
            {
                textBox.Background = Brushes.White;
            }
            else
            {
                textBox.Background = Brushes.Red;
            }
        }

        private void INT_TextBox_TextChanged(object sender, TextChangedEventArgs e)
        {
            var textBox = sender as TextBox;
            int var;

            if (int.TryParse(textBox.Text, NumberStyles.Number, CultureInfo.InvariantCulture, out var))
            {
                textBox.Background = Brushes.White;
            }
            else
            {
                textBox.Background = Brushes.Red;
            }
        }

        private void Generate_click_button(object sender, RoutedEventArgs e)
        {
            this.DialogResult = true;
        }
    }
}
