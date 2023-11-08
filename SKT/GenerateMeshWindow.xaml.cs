using System;
using System.Collections.Generic;
using System.Collections.ObjectModel;
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
    /// Логика взаимодействия для GenerateMeshWindow.xaml
    /// </summary>
    public partial class GenerateMeshWindow : Window
    {
        public GenerateMeshWindow()
        {
            InitializeComponent();
            DataContext = this;
        }

        int nx = 3, nz = 3;
        double min_x = -200, max_x = 200, min_z = -200, max_z = -100;

        public int NX { get { return nx; } set { if (value > 1) nx = value; } }
        public int NZ { get { return nz; } set { if (value > 1) nz = value; } }

        public double MinX { get { return min_x; } set { if (value < max_x) min_x = value; } }

        public double MaxX { get { return max_x; } set { if (value > min_x) max_x = value; } }

        public double MinZ { get { return min_z; } set { if (value < max_z) min_z = value; } }

        public double MaxZ { get { return max_z; } set { if (value > min_z) max_z = value; } }


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
