using System;
using System.CodeDom;
using System.Collections.Generic;
using System.ComponentModel;
using System.Globalization;
using System.Linq;
using System.Windows;
using System.Windows.Media.Media3D;
using static SKT.Solver;

namespace SKT
{
    public enum CalculationMode
    {
        Center,
        Gauss2,
        Gauss3,
        Gauss7,
    }

    public class InverseInput
    {
        double start_gamma;
        double allow_max_difference;
        double gamma_mult;
        double max_functional_in_percent;
        int max_iterations;

        public double StartGamma
        {
            get { return start_gamma; }
            set
            {
                if (value > 0)
                {
                    start_gamma = value;
                }
            }
        }

        public double AllowMaxDifference
        {
            get { return allow_max_difference; }
            set
            {
                if (value > 1.0)
                {
                    allow_max_difference = value;
                }
            }
        }

        public double GammaMult
        {
            get { return gamma_mult; }
            set
            {
                if (value > 1.0)
                {
                    gamma_mult = value;
                }
            }
        }

        public double MaxFunctionalInPercent
        {
            get { return max_functional_in_percent; }
            set
            {
                if (value > 0)
                {
                    max_functional_in_percent = value;
                }
            }
        }

        public int MaxIterations
        {
            get { return max_iterations; }
            set
            {
                if (value > 0)
                {
                    max_iterations = value;
                }
            }
        }

        public bool SaveAllIterations { get; set; }

        public InverseInput()
        {
            StartGamma = 1.0e-9;
            AllowMaxDifference = 2.0;
            GammaMult = 10.0;
            MaxFunctionalInPercent = 200.0;
            MaxIterations = 10;
            SaveAllIterations = true;
        }
    }

    public class InverseIterationResult
    {
        // Намагниченность в ячейках после выполнения инверсии
        public List<double> ValuesP { get; set; }

        // Значение поля в точках измерения после инверсии
        public List<double[]> ValuesB { get; set; }

        // Значение функционала после инверсий
        public double Func { get; set; }

        // Номера итерацим
        public int NumIter { get; set; }

        public InverseIterationResult()
        {
            ValuesP = new();
            ValuesB = new();
            Func = 0.0;
            NumIter = 0;
        }
    }

    public class InverseResult
    {
        public List<InverseIterationResult> Iters { get; set; }

        public InverseResult()
        {
            Iters = new();
        }

        public void Clear()
        {
            Iters.Clear();
        }
    }

    public class Solver
    {
        // на вход подаются сетка и профили
        // на выходе заполняется значение поля в точках измерения в профилях
        public static void SolveForward(
            RectangleMesh mesh,
            List<Profile_receivers> profiles,
            CalculationMode mode,
            double[]? field_p_in_mesh_elements = null,
            List<double[]>? field_B_in_receivers = null)
        {
            // Глобальные координаты точки Гаусса или центра
            Point3D glob_coord = new Point3D();

            // Координаты точки измерения, при переносе начала координат в точку Гаусса или центр
            Point new_coord = new();

            // Площадь элемента
            double mes = mesh.Mes;

            // Расстояние от точки измерений до k-ой точки Гаусса или центра
            double rk_1, rk_2, rk_3;

            double I = 1.0;

            double value_k;
            double value;

            switch (mode)
            {
                case CalculationMode.Center:
                    for (int i_profile = 0; i_profile < profiles.Count; i_profile++)
                    {
                        var profile = profiles[i_profile];
                        for (int i = 0; i < profile.Coords.Count; i++)
                        {
                            Point old_coord_receiver = profile.Coords[i];
                            value = 0.0;

                            for (int i_x = 0; i_x < mesh.X.Count - 1; i_x++)
                            {
                                for (int j_z = 0; j_z < mesh.Z.Count - 1; j_z++)
                                {
                                    double x0 = mesh.X[i_x];
                                    double x1 = mesh.X[i_x + 1];

                                    double z0 = mesh.Z[j_z];
                                    double z1 = mesh.Z[j_z + 1];

                                    double P;
                                    if (field_p_in_mesh_elements != null)
                                    {
                                        P = field_p_in_mesh_elements[i_x * (mesh.Z.Count - 1) + j_z];
                                    }
                                    else
                                    {
                                        P = mesh.P[i_x * (mesh.Z.Count - 1) + j_z].Value;
                                    }

                                    glob_coord.X = 0.5 * (x1 - x0) + x0;
                                    glob_coord.Z = 0.5 * (z1 - z0) + z0;

                                    new_coord.X = old_coord_receiver.X - glob_coord.X;
                                    new_coord.Z = old_coord_receiver.Z - glob_coord.Z;

                                    rk_1 = Math.Sqrt(new_coord.X * new_coord.X + new_coord.Z * new_coord.Z);
                                    rk_2 = rk_1 * rk_1;
                                    rk_3 = rk_2 * rk_1;

                                    value += mes * I / (4.0 * Math.PI * rk_3) * (P * (3.0 * new_coord.X * new_coord.X / rk_2 - 1.0));
                                }
                            }

                            if (field_B_in_receivers != null)
                            {
                                field_B_in_receivers[i_profile][i] = value;
                            }
                            else
                            {
                                profile.Values[i] = value;
                            }
                        }
                    }
                    break;

                case CalculationMode.Gauss2:
                    // Квадратура Гаусса-2
                    double gauss2_ksi_1 = -1.0 / Math.Sqrt(3);
                    double gauss2_ksi_2 = 1.0 / Math.Sqrt(3);

                    double[] gauss2_loc_x = { gauss2_ksi_1, gauss2_ksi_2, gauss2_ksi_1, gauss2_ksi_2 };
                    double[] gauss2_loc_y = { gauss2_ksi_1, gauss2_ksi_1, gauss2_ksi_2, gauss2_ksi_2 };
                    double[] gauss2_w = { 1.0, 1.0, 1.0, 1.0, };

                    for (int i_profile = 0; i_profile < profiles.Count; i_profile++)
                    {
                        var profile = profiles[i_profile];

                        for (int i = 0; i < profile.Coords.Count; i++)
                        {
                            Point old_coord_receiver = profile.Coords[i];
                            value = 0.0;

                            for (int i_x = 0; i_x < mesh.X.Count - 1; i_x++)
                            {
                                for (int j_z = 0; j_z < mesh.Z.Count - 1; j_z++)
                                {
                                    double x0 = mesh.X[i_x];
                                    double x1 = mesh.X[i_x + 1];

                                    double z0 = mesh.Z[j_z];
                                    double z1 = mesh.Z[j_z + 1];

                                    double P;
                                    if (field_p_in_mesh_elements != null)
                                    {
                                        P = field_p_in_mesh_elements[i_x * (mesh.Z.Count - 1) + j_z];
                                    }
                                    else
                                    {
                                        P = mesh.P[i_x * (mesh.Z.Count - 1) + j_z].Value;
                                    }

                                    value_k = 0.0;

                                    for (int k = 0; k < 4; k++)
                                    {
                                        glob_coord.X = 0.5 * (x1 - x0) * (gauss2_loc_x[k] + 1.0) + x0;
                                        glob_coord.Z = 0.5 * (z1 - z0) * (gauss2_loc_y[k] + 1.0) + z0;

                                        new_coord.X = old_coord_receiver.X - glob_coord.X;
                                        new_coord.Z = old_coord_receiver.Z - glob_coord.Z;

                                        rk_1 = Math.Sqrt(new_coord.X * new_coord.X + new_coord.Z * new_coord.Z);
                                        rk_2 = rk_1 * rk_1;
                                        rk_3 = rk_2 * rk_1;

                                        value_k += gauss2_w[k] / (4.0 * Math.PI * rk_3) * (P * (3.0 * new_coord.X * new_coord.X / rk_2 - 1.0));
                                    }

                                    value += value_k * mes * I;
                                }
                            }

                            if (field_B_in_receivers != null)
                            {
                                field_B_in_receivers[i_profile][i] = value;
                            }
                            else
                            {
                                profile.Values[i] = value;
                            }
                        }
                    }
                    break;

                case CalculationMode.Gauss3:
                    // Квадратура Гаусса-3
                    double gauss3_ksi_1 = -Math.Sqrt(3.0 / 5.0);
                    double gauss3_ksi_2 = 0;
                    double gauss3_ksi_3 = Math.Sqrt(3.0 / 5.0);

                    double gauss3_w1 = 5.0 / 9.0;
                    double gauss3_w2 = 8.0 / 9.0;
                    double gauss3_w3 = 5.0 / 9.0;

                    double[] gauss3_loc_x = { 
                        gauss3_ksi_1, gauss3_ksi_2, gauss3_ksi_3,
                        gauss3_ksi_1, gauss3_ksi_2, gauss3_ksi_3,
                        gauss3_ksi_1, gauss3_ksi_2, gauss3_ksi_3
                    };
                    double[] gauss3_loc_y = {
                        gauss3_ksi_1, gauss3_ksi_1, gauss3_ksi_1,
                        gauss3_ksi_2, gauss3_ksi_2, gauss3_ksi_2,
                        gauss3_ksi_3, gauss3_ksi_3, gauss3_ksi_3
                    };
                    double[] gauss3_w = {
                        gauss3_w1 * gauss3_w1, gauss3_w1 * gauss3_w2, gauss3_w1 * gauss3_w3,
                        gauss3_w2 * gauss3_w1, gauss3_w2 * gauss3_w2, gauss3_w2 * gauss3_w3,
                        gauss3_w3 * gauss3_w1, gauss3_w3 * gauss3_w2, gauss3_w3 * gauss3_w3
                    };

                    for (int i_profile = 0; i_profile < profiles.Count; i_profile++)
                    {
                        var profile = profiles[i_profile];

                        for (int i = 0; i < profile.Coords.Count; i++)
                        {
                            Point old_coord_receiver = profile.Coords[i];
                            value = 0.0;

                            for (int i_x = 0; i_x < mesh.X.Count - 1; i_x++)
                            {
                                for (int j_z = 0; j_z < mesh.Z.Count - 1; j_z++)
                                {
                                    double x0 = mesh.X[i_x];
                                    double x1 = mesh.X[i_x + 1];

                                    double z0 = mesh.Z[j_z];
                                    double z1 = mesh.Z[j_z + 1];

                                    double P;
                                    if (field_p_in_mesh_elements != null)
                                    {
                                        P = field_p_in_mesh_elements[i_x * (mesh.Z.Count - 1) + j_z];
                                    }
                                    else
                                    {
                                        P = mesh.P[i_x * (mesh.Z.Count - 1) + j_z].Value;
                                    }

                                    value_k = 0.0;

                                    for (int k = 0; k < 9; k++)
                                    {
                                        glob_coord.X = 0.5 * (x1 - x0) * (gauss3_loc_x[k] + 1.0) + x0;
                                        glob_coord.Z = 0.5 * (z1 - z0) * (gauss3_loc_y[k] + 1.0) + z0;

                                        new_coord.X = old_coord_receiver.X - glob_coord.X;
                                        new_coord.Z = old_coord_receiver.Z - glob_coord.Z;

                                        rk_1 = Math.Sqrt(new_coord.X * new_coord.X + new_coord.Z * new_coord.Z);
                                        rk_2 = rk_1 * rk_1;
                                        rk_3 = rk_2 * rk_1;

                                        value_k += gauss3_w[k] / (4.0 * Math.PI * rk_3) * (P * (3.0 * new_coord.X * new_coord.X / rk_2 - 1.0));
                                    }

                                    value += value_k * mes * I;
                                }
                            }

                            if (field_B_in_receivers != null)
                            {
                                field_B_in_receivers[i_profile][i] = value;
                            }
                            else
                            {
                                profile.Values[i] = value;
                            }
                        }
                    }
                    break;

                case CalculationMode.Gauss7:
                    // Квадратура Гаусса для четырехугольного мастер-элемента 7-го порядка аппроксимации
                    double a = Math.Sqrt((114.0 - 3.0 * Math.Sqrt(583.0)) / 287.0);
                    double b = Math.Sqrt((114.0 + 3.0 * Math.Sqrt(583.0)) / 287.0);
                    double c = Math.Sqrt(6.0 / 7.0);
                    double wa = 307.0 / 810.0 + 923.0 / (270.0 * Math.Sqrt(583.0));
                    double wb = 307.0 / 810.0 - 923.0 / (270.0 * Math.Sqrt(583.0));
                    double wc = 98.0 / 405.0;

                    double[] loc_x = { -c, c, 0.0, 0.0, -a, a, -a, a, -b, b, -b, b };
                    double[] loc_y = { 0.0, 0.0, -c, c, -a, -a, a, a, -b, -b, b, b };
                    double[] w = { wa, wa, wa, wa, wc, wc, wc, wc, wb, wb, wb, wb };

                    for (int i_profile = 0; i_profile < profiles.Count; i_profile++)
                    {
                        var profile = profiles[i_profile];

                        for (int i = 0; i < profile.Coords.Count; i++)
                        {
                            Point old_coord_receiver = profile.Coords[i];
                            value = 0.0;

                            for (int i_x = 0; i_x < mesh.X.Count - 1; i_x++)
                            {
                                for (int j_z = 0; j_z < mesh.Z.Count - 1; j_z++)
                                {
                                    double x0 = mesh.X[i_x];
                                    double x1 = mesh.X[i_x + 1];

                                    double z0 = mesh.Z[j_z];
                                    double z1 = mesh.Z[j_z + 1];

                                    double P;
                                    if (field_p_in_mesh_elements != null)
                                    {
                                        P = field_p_in_mesh_elements[i_x * (mesh.Z.Count - 1) + j_z];
                                    }
                                    else
                                    {
                                        P = mesh.P[i_x * (mesh.Z.Count - 1) + j_z].Value;
                                    }

                                    value_k = 0.0;

                                    for (int k = 0; k < 12; k++)
                                    {
                                        glob_coord.X = 0.5 * (x1 - x0) * (loc_x[k] + 1.0) + x0;
                                        glob_coord.Z = 0.5 * (z1 - z0) * (loc_y[k] + 1.0) + z0;

                                        new_coord.X = old_coord_receiver.X - glob_coord.X;
                                        new_coord.Z = old_coord_receiver.Z - glob_coord.Z;

                                        rk_1 = Math.Sqrt(new_coord.X * new_coord.X + new_coord.Z * new_coord.Z);
                                        rk_2 = rk_1 * rk_1;
                                        rk_3 = rk_2 * rk_1;

                                        value_k += w[k] / (4.0 * Math.PI * rk_3) * (P * (3.0 * new_coord.X * new_coord.X / rk_2 - 1.0));
                                    }

                                    value += value_k * mes * I;
                                }
                            }

                            if (field_B_in_receivers != null)
                            {
                                field_B_in_receivers[i_profile][i] = value;
                            }
                            else
                            {
                                profile.Values[i] = value;
                            }
                        }
                    }
                    break;

                default:
                    break;
            }
        }






        public static InverseResult SolveInverse(
            RectangleMesh mesh, 
            List<Profile_receivers> profiles,
            CalculationMode mode,
            double start_gamma,
            double allow_max_difference, 
            double gamma_mult,
            double max_functional_in_percent,
            int max_iterations,
            bool save_all_iterations)
        {
            InverseResult result = new();

            int i, j;

            double[][] gamma = MatrixCreate(mesh.X.Count - 1, mesh.Z.Count - 1);
            for (i = 0; i < mesh.X.Count - 1; i++)
            {
                for (j = 0; j < mesh.Z.Count - 1; j++)
                {
                    gamma[i][j] = start_gamma;
                }
            }

            double first_func = 0.0;
            double current_func;

            bool is_first_inverse = true;
            bool is_end = false;

            int iter = 0;

            do
            {
                iter++;

                double[] field_p_in_mesh_elements = new double[mesh.CountElements];

                List<double[]> field_B_in_profiles = new List<double[]>();
                for (i = 0; i < profiles.Count; i++)
                {
                    double[] field_B_in_profile = new double[profiles[i].Count];
                    field_B_in_profiles.Add(field_B_in_profile);
                }

                Solver.SolveInverse(mesh, profiles, mode, field_p_in_mesh_elements, gamma);
                Solver.SolveForward(mesh, profiles, mode, field_p_in_mesh_elements, field_B_in_profiles);

                current_func = CalculateFunctinal(mesh, profiles, gamma, field_p_in_mesh_elements, field_B_in_profiles);


                if (is_first_inverse)
                {
                    first_func = current_func;

                    InverseIterationResult iter_result = new();

                    iter_result.ValuesP = field_p_in_mesh_elements.ToList();
                    iter_result.ValuesB = field_B_in_profiles;
                    iter_result.Func = current_func;
                    iter_result.NumIter = iter;

                    result.Iters.Add(iter_result);
                }

                if (((current_func / first_func - 1.0) * 100.0 >= max_functional_in_percent) || (iter >= max_iterations))
                {
                    is_end = true;

                    if (!is_first_inverse)
                    {
                        InverseIterationResult iter_result = new();

                        iter_result.ValuesP = field_p_in_mesh_elements.ToList();
                        iter_result.ValuesB = field_B_in_profiles;
                        iter_result.Func = current_func;
                        iter_result.NumIter = iter;

                        result.Iters.Add(iter_result);
                    }
                }

                if (!is_end)
                {
                    UpdateGamma(gamma, mesh.X.Count - 1, mesh.Z.Count - 1, field_p_in_mesh_elements, gamma_mult, allow_max_difference);
                }

                if (save_all_iterations && !is_end && !is_first_inverse)
                {
                    InverseIterationResult iter_result = new();

                    iter_result.ValuesP = field_p_in_mesh_elements.ToList();
                    iter_result.ValuesB = field_B_in_profiles;
                    iter_result.Func = current_func;
                    iter_result.NumIter = iter;

                    result.Iters.Add(iter_result);
                }

                is_first_inverse = false;

            } while (!is_end);

            return result;
        }







        // Если filed_p_in_mesh_elements == null, то найденные значения записаны в mesh.P, иначе в filed_p_in_mesh_elements
        // gamma != null, то будет выполнена регуляризация
        public static void SolveInverse(
            RectangleMesh mesh,
            List<Profile_receivers> receivers,
            CalculationMode mode,
            double[]? field_p_in_mesh_elements = null,
            double[][]? gamma = null)
        {
            // создать матрицу L, A
            // заполнить L в зависимости от выбранных точек гаусса
            // перемножить и получить A
            // добавить регуляризацию
            // решить слау

            // K - ячеек, n - примеников
            // n = количество профилей * количество точек в каждом прифиле
            // K = количетсво ячеек
            // L - n строк, K столбцов

            // S - вектор исходных данных, длина равна n
            // P - вектор параметров, которые нужно найти длина - K

            // A = LT * L        [K, n] * [n, K] = [K, K]
            // b = LT * S        [K, n] * [n] = [K]
            // Ap = b            [K, K] * [K] = [K]

            // n
            int rows_l = receivers.Count * receivers[0].Count;

            // K
            int cols_l = mesh.CountElements;

            double[][] L = MatrixCreate(rows_l, cols_l);
            double[][] A = MatrixCreate(cols_l, cols_l);
            double[] b = new double[cols_l];
            double[] S = new double[rows_l];
            double[] p = new double[cols_l];
            double[] solution = new double[cols_l];

            double mes = mesh.Mes;

            double new_coord_x;
            double new_coord_z;

            double rk_1;
            double rk_2;
            double rk_3;

            double x0;
            double x1;

            double z0;
            double z1;

            switch (mode)
            {
                case CalculationMode.Center:
                    for (int index_receivers = 0, index_n = 0; index_receivers < receivers.Count; index_receivers++)
                    {
                        for (int index_coords = 0; index_coords < receivers[index_receivers].Coords.Count; index_coords++, index_n++)
                        {
                            Point old_coord_receiver = receivers[index_receivers].Coords[index_coords];

                            S[index_n] = receivers[index_receivers].Values[index_coords];

                            for (int i_x = 0, index_K = 0; i_x < mesh.X.Count - 1; i_x++)
                            {
                                for (int j_z = 0; j_z < mesh.Z.Count - 1; j_z++, index_K++)
                                {
                                    x0 = mesh.X[i_x];
                                    x1 = mesh.X[i_x + 1];

                                    z0 = mesh.Z[j_z];
                                    z1 = mesh.Z[j_z + 1];

                                    new_coord_x = old_coord_receiver.X - (0.5 * (x1 - x0) + x0);
                                    new_coord_z = old_coord_receiver.Z - (0.5 * (z1 - z0) + z0);

                                    rk_1 = Math.Sqrt(new_coord_x * new_coord_x + new_coord_z * new_coord_z);
                                    rk_2 = rk_1 * rk_1;
                                    rk_3 = rk_2 * rk_1;

                                    L[index_n][index_K] = mes / (4.0 * Math.PI * rk_3) * (3.0 * new_coord_x * new_coord_x / rk_2 - 1.0);
                                }
                            }
                        }
                    }
                    break;

                case CalculationMode.Gauss2:
                    double gauss2_ksi_1 = -1.0 / Math.Sqrt(3);
                    double gauss2_ksi_2 = 1.0 / Math.Sqrt(3);

                    double[] gauss2_loc_x = { gauss2_ksi_1, gauss2_ksi_2, gauss2_ksi_1, gauss2_ksi_2 };
                    double[] gauss2_loc_y = { gauss2_ksi_1, gauss2_ksi_1, gauss2_ksi_2, gauss2_ksi_2 };
                    double[] gauss2_w = { 1.0, 1.0, 1.0, 1.0, };

                    double value;

                    for (int index_receivers = 0, index_n = 0; index_receivers < receivers.Count; index_receivers++)
                    {
                        for (int index_coords = 0; index_coords < receivers[index_receivers].Coords.Count; index_coords++, index_n++)
                        {
                            Point old_coord_receiver = receivers[index_receivers].Coords[index_coords];

                            S[index_n] = receivers[index_receivers].Values[index_coords];

                            for (int i_x = 0, index_K = 0; i_x < mesh.X.Count - 1; i_x++)
                            {
                                for (int j_z = 0; j_z < mesh.Z.Count - 1; j_z++, index_K++)
                                {
                                    x0 = mesh.X[i_x];
                                    x1 = mesh.X[i_x + 1];

                                    z0 = mesh.Z[j_z];
                                    z1 = mesh.Z[j_z + 1];

                                    value = 0.0;

                                    for (int k = 0; k < 4; k++)
                                    {
                                        new_coord_x = old_coord_receiver.X - (0.5 * (x1 - x0) * (gauss2_loc_x[k] + 1.0) + x0);
                                        new_coord_z = old_coord_receiver.Z - (0.5 * (z1 - z0) * (gauss2_loc_y[k] + 1.0) + z0);

                                        rk_1 = Math.Sqrt(new_coord_x * new_coord_x + new_coord_z * new_coord_z);
                                        rk_2 = rk_1 * rk_1;
                                        rk_3 = rk_2 * rk_1;

                                        value += gauss2_w[k] / (4.0 * Math.PI * rk_3) * (3.0 * new_coord_x * new_coord_x / rk_2 - 1.0);
                                    }

                                    L[index_n][index_K] = mes * value;
                                }
                            }
                        }
                    }
                    break;

                case CalculationMode.Gauss3:
                    double gauss3_ksi_1 = -Math.Sqrt(3.0 / 5.0);
                    double gauss3_ksi_2 = 0;
                    double gauss3_ksi_3 = Math.Sqrt(3.0 / 5.0);

                    double gauss3_w1 = 5.0 / 9.0;
                    double gauss3_w2 = 8.0 / 9.0;
                    double gauss3_w3 = 5.0 / 9.0;

                    double[] gauss3_loc_x = {
                        gauss3_ksi_1, gauss3_ksi_2, gauss3_ksi_3,
                        gauss3_ksi_1, gauss3_ksi_2, gauss3_ksi_3,
                        gauss3_ksi_1, gauss3_ksi_2, gauss3_ksi_3
                    };
                    double[] gauss3_loc_y = {
                        gauss3_ksi_1, gauss3_ksi_1, gauss3_ksi_1,
                        gauss3_ksi_2, gauss3_ksi_2, gauss3_ksi_2,
                        gauss3_ksi_3, gauss3_ksi_3, gauss3_ksi_3
                    };
                    double[] gauss3_w = {
                        gauss3_w1 * gauss3_w1, gauss3_w1 * gauss3_w2, gauss3_w1 * gauss3_w3,
                        gauss3_w2 * gauss3_w1, gauss3_w2 * gauss3_w2, gauss3_w2 * gauss3_w3,
                        gauss3_w3 * gauss3_w1, gauss3_w3 * gauss3_w2, gauss3_w3 * gauss3_w3
                    };

                    for (int index_receivers = 0, index_n = 0; index_receivers < receivers.Count; index_receivers++)
                    {
                        for (int index_coords = 0; index_coords < receivers[index_receivers].Coords.Count; index_coords++, index_n++)
                        {
                            Point old_coord_receiver = receivers[index_receivers].Coords[index_coords];

                            S[index_n] = receivers[index_receivers].Values[index_coords];

                            for (int i_x = 0, index_K = 0; i_x < mesh.X.Count - 1; i_x++)
                            {
                                for (int j_z = 0; j_z < mesh.Z.Count - 1; j_z++, index_K++)
                                {
                                    x0 = mesh.X[i_x];
                                    x1 = mesh.X[i_x + 1];

                                    z0 = mesh.Z[j_z];
                                    z1 = mesh.Z[j_z + 1];

                                    value = 0.0;

                                    for (int k = 0; k < 9; k++)
                                    {
                                        new_coord_x = old_coord_receiver.X - (0.5 * (x1 - x0) * (gauss3_loc_x[k] + 1.0) + x0);
                                        new_coord_z = old_coord_receiver.Z - (0.5 * (z1 - z0) * (gauss3_loc_y[k] + 1.0) + z0);

                                        rk_1 = Math.Sqrt(new_coord_x * new_coord_x + new_coord_z * new_coord_z);
                                        rk_2 = rk_1 * rk_1;
                                        rk_3 = rk_2 * rk_1;

                                        value += gauss3_w[k] / (4.0 * Math.PI * rk_3) * (3.0 * new_coord_x * new_coord_x / rk_2 - 1.0);
                                    }

                                    L[index_n][index_K] = mes * value;
                                }
                            }
                        }
                    }
                    break;

                case CalculationMode.Gauss7:
                    double gauss_a = Math.Sqrt((114.0 - 3.0 * Math.Sqrt(583.0)) / 287.0);
                    double gauss_b = Math.Sqrt((114.0 + 3.0 * Math.Sqrt(583.0)) / 287.0);
                    double gauss_c = Math.Sqrt(6.0 / 7.0);
                    double gauss_wa = 307.0 / 810.0 + 923.0 / (270.0 * Math.Sqrt(583.0));
                    double gauss_wb = 307.0 / 810.0 - 923.0 / (270.0 * Math.Sqrt(583.0));
                    double gauss_wc = 98.0 / 405.0;

                    double[] gauss7_loc_x = { -gauss_c, gauss_c, 0.0, 0.0, -gauss_a, gauss_a, -gauss_a, gauss_a, -gauss_b, gauss_b, -gauss_b, gauss_b };
                    double[] gauss7_loc_y = { 0.0, 0.0, -gauss_c, gauss_c, -gauss_a, -gauss_a, gauss_a, gauss_a, -gauss_b, -gauss_b, gauss_b, gauss_b };
                    double[] gauss7_w = { gauss_wa, gauss_wa, gauss_wa, gauss_wa, gauss_wc, gauss_wc, gauss_wc, gauss_wc, gauss_wb, gauss_wb, gauss_wb, gauss_wb };

                    for (int index_receivers = 0, index_n = 0; index_receivers < receivers.Count; index_receivers++)
                    {
                        for (int index_coords = 0; index_coords < receivers[index_receivers].Coords.Count; index_coords++, index_n++)
                        {
                            Point old_coord_receiver = receivers[index_receivers].Coords[index_coords];

                            S[index_n] = receivers[index_receivers].Values[index_coords];

                            for (int i_x = 0, index_K = 0; i_x < mesh.X.Count - 1; i_x++)
                            {
                                for (int j_z = 0; j_z < mesh.Z.Count - 1; j_z++, index_K++)
                                {
                                    x0 = mesh.X[i_x];
                                    x1 = mesh.X[i_x + 1];

                                    z0 = mesh.Z[j_z];
                                    z1 = mesh.Z[j_z + 1];

                                    value = 0.0;

                                    for (int k = 0; k < 12; k++)
                                    {
                                        new_coord_x = old_coord_receiver.X - (0.5 * (x1 - x0) * (gauss7_loc_x[k] + 1.0) + x0);
                                        new_coord_z = old_coord_receiver.Z - (0.5 * (z1 - z0) * (gauss7_loc_y[k] + 1.0) + z0);

                                        rk_1 = Math.Sqrt(new_coord_x * new_coord_x + new_coord_z * new_coord_z);
                                        rk_2 = rk_1 * rk_1;
                                        rk_3 = rk_2 * rk_1;

                                        value += gauss7_w[k] / (4.0 * Math.PI * rk_3) * (3.0 * new_coord_x * new_coord_x / rk_2 - 1.0);
                                    }

                                    L[index_n][index_K] = mes * value;
                                }
                            }
                        }
                    }
                    break;

                default:
                    break;
            }

            // A = LT * L        [K, n] * [n, K] = [K, K]
            for (int k = 0; k < rows_l; k++)
            {
                double[] temp = L[k];

                for (int i = 0; i < cols_l; i++)
                {
                    double[] temp_2 = A[i];
                    double temp_3 = temp[i];

                    for (int j = 0; j < cols_l; j++)
                    {
                        temp_2[j] += temp_3 * temp[j];
                    }
                }
            }

            // b = LT * S        [K, n] * [n] = [K]
            for (int j = 0; j < rows_l; j++)
            {
                double[] temp = L[j];
                double temp_1 = S[j];

                for (int i = 0; i < cols_l; i++)
                {
                    b[i] += temp[i] * temp_1;
                }
            }

            // РЕГУЛЯРИЗАЦИЯ

            if (gamma != null)
            {
                int current_index, neighbor_index;
                int size_z = mesh.Z.Count - 1;
                int size_x = mesh.X.Count - 1;
                double gamma_ij;

                for (int i = 0; i < mesh.X.Count - 1; i++)
                {
                    for (int j = 0; j < mesh.Z.Count - 1; j++)
                    {
                        current_index = i * size_z + j;
                        gamma_ij = gamma[i][j];

                        if (i != 0)
                        {
                            if (j != 0)
                            {
                                neighbor_index = (i - 1) * size_z + j - 1;
                                A[current_index][neighbor_index] -= (gamma_ij + gamma[i - 1][j - 1]);
                                A[current_index][current_index] += gamma_ij + gamma[i - 1][j - 1];
                            }

                            neighbor_index = (i - 1) * size_z + j;
                            A[current_index][neighbor_index] -= (gamma_ij + gamma[i - 1][j]);
                            A[current_index][current_index] += gamma_ij + gamma[i - 1][j];

                            if (j != size_z - 1)
                            {
                                neighbor_index = (i - 1) * size_z + j + 1;
                                A[current_index][neighbor_index] -= (gamma_ij + gamma[i - 1][j + 1]);
                                A[current_index][current_index] += gamma_ij + gamma[i - 1][j + 1];
                            }
                        }

                        if (j != 0)
                        {
                            neighbor_index = i * size_z + j - 1;
                            A[current_index][neighbor_index] -= (gamma_ij + gamma[i][j - 1]);
                            A[current_index][current_index] += gamma_ij + gamma[i][j - 1];
                        }

                        if (j != size_z - 1)
                        {
                            neighbor_index = i * size_z + j + 1;
                            A[current_index][neighbor_index] -= (gamma_ij + gamma[i][j + 1]);
                            A[current_index][current_index] += gamma_ij + gamma[i][j + 1];
                        }

                        if (i != size_x - 1)
                        {
                            if (j != 0)
                            {
                                neighbor_index = (i + 1) * size_z + j - 1;
                                A[current_index][neighbor_index] -= (gamma_ij + gamma[i + 1][j - 1]);
                                A[current_index][current_index] += gamma_ij + gamma[i + 1][j - 1];
                            }

                            neighbor_index = (i + 1) * size_z + j;
                            A[current_index][neighbor_index] -= (gamma_ij + gamma[i + 1][j]);
                            A[current_index][current_index] += gamma_ij + gamma[i + 1][j];

                            if (j != size_z - 1)
                            {
                                neighbor_index = (i + 1) * size_z + j + 1;
                                A[current_index][neighbor_index] -= (gamma_ij + gamma[i + 1][j + 1]);
                                A[current_index][current_index] += gamma_ij + gamma[i + 1][j + 1];
                            }
                        }
                    }
                }
            }

            // Решение СЛАУ Ap=b

            Gauss(A, b, cols_l);

            if (field_p_in_mesh_elements != null)
            {
                for (int i_x = 0, index_K = 0; i_x < mesh.X.Count - 1; i_x++)
                {
                    for (int j_z = 0; j_z < mesh.Z.Count - 1; j_z++, index_K++)
                    {
                        field_p_in_mesh_elements[index_K] = b[index_K];
                    }
                }
            }
            else
            {
                for (int i_x = 0, index_K = 0; i_x < mesh.X.Count - 1; i_x++)
                {
                    for (int j_z = 0; j_z < mesh.Z.Count - 1; j_z++, index_K++)
                    {
                        mesh.P[index_K].Value = b[index_K];
                    }
                }
            }
        }






        public static double[][] MatrixCreate(int rows, int cols)
        {
            double[][] result = new double[rows][];
            for (int i = 0; i < rows; ++i)
            {
                result[i] = new double[cols];
            }
            return result;
        }






        public static void Gauss(double[][] A, double[] b, int N)
        {
            int i, j, k, n, m;
            n = N;
            double aa, bb;
            for (k = 0; k < n; k++) //Поиск максимального элемента в первом столбце
            {
                aa = Math.Abs(A[k][k]);
                i = k;
                for (m = k + 1; m < n; m++)
                {
                    if (Math.Abs(A[m][k]) > aa)
                    {
                        i = m;
                        aa = Math.Abs(A[m][k]);
                    }
                }

                if (aa == 0)   //проверка на нулевой элемент
                {
                    MessageBox.Show("Система не имеет решений");
                }

                if (i != k)  //  перестановка i-ой строки, содержащей главный элемент k-ой строки
                {
                    for (j = k; j < n; j++)
                    {
                        bb = A[k][j];
                        A[k][j] = A[i][j];
                        A[i][j] = bb;
                    }
                    bb = b[k];
                    b[k] = b[i];
                    b[i] = bb;
                }
                aa = A[k][k];//преобразование k-ой строки (Вычисление масштабирующих множителей)
                A[k][k] = 1;
                for (j = k + 1; j < n; j++)
                {
                    A[k][j] = A[k][j] / aa;
                }
                b[k] /= aa;

                //преобразование строк с помощью k-ой строки
                for (i = k + 1; i < n; i++)
                {
                    bb = A[i][k];
                    A[i][k] = 0;

                    if (bb != 0)
                    {
                        for (j = k + 1; j < n; j++)
                        {
                            A[i][j] = A[i][j] - bb * A[k][j];
                        }

                        b[i] -= bb * b[k];
                    }

                }
            }

            for (i = n - 1; i >= 0; i--)   //Нахождение решений СЛАУ
            {
                for (j = n - 1; j > i; j--)
                    b[i] -= A[i][j] * b[j];
            }
        }






        public static double CalculateFunctinal(
            RectangleMesh mesh,
            List<Profile_receivers> profiles,
            double[][] gamma,
            double[] field_p_in_mesh_elements,
            List<double[]> field_B_in_profiles)
        {
            double functional = 0.0;

            int i, j;
            int current_index;
            int neighbor_index;
            double difference;
            double gamma_ij;

            for (int i_profile = 0; i_profile < profiles.Count; i_profile++)
            {
                var profile = profiles[i_profile];
                var field_B_in_profile = field_B_in_profiles[i_profile];
                for (i = 0; i < profile.Coords.Count; i++)
                {
                    functional += (profile.Values[i] - field_B_in_profile[i]) * (profile.Values[i] - field_B_in_profile[i]);
                }
            }

            // внутренние элементы
            for (i = 1; i < mesh.X.Count - 2; i++)
            {
                for (j = 1; j < mesh.Z.Count - 2; j++)
                {
                    current_index = i * (mesh.Z.Count - 1) + j;
                    gamma_ij = gamma[i][j];

                    neighbor_index = (i - 1) * (mesh.Z.Count - 1) + j - 1;
                    difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                    functional += gamma_ij * difference * difference;

                    neighbor_index = (i - 1) * (mesh.Z.Count - 1) + j;
                    difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                    functional += gamma_ij * difference * difference;

                    neighbor_index = (i - 1) * (mesh.Z.Count - 1) + j + 1;
                    difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                    functional += gamma_ij * difference * difference;

                    neighbor_index = i * (mesh.Z.Count - 1) + j - 1;
                    difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                    functional += gamma_ij * difference * difference;

                    neighbor_index = i * (mesh.Z.Count - 1) + j + 1;
                    difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                    functional += gamma_ij * difference * difference;

                    neighbor_index = (i + 1) * (mesh.Z.Count - 1) + j - 1;
                    difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                    functional += gamma_ij * difference * difference;

                    neighbor_index = (i + 1) * (mesh.Z.Count - 1) + j;
                    difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                    functional += gamma_ij * difference * difference;

                    neighbor_index = (i + 1) * (mesh.Z.Count - 1) + j + 1;
                    difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                    functional += gamma_ij * difference * difference;
                }
            }

            // элементы на нижней границе
            j = 0;
            for (i = 1; i < mesh.X.Count - 2; i++)
            {
                current_index = i * (mesh.Z.Count - 1) + j;
                gamma_ij = gamma[i][j];

                neighbor_index = (i - 1) * (mesh.Z.Count - 1) + j;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;

                neighbor_index = (i - 1) * (mesh.Z.Count - 1) + j + 1;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;

                neighbor_index = i * (mesh.Z.Count - 1) + j + 1;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;

                neighbor_index = (i + 1) * (mesh.Z.Count - 1) + j;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;

                neighbor_index = (i + 1) * (mesh.Z.Count - 1) + j + 1;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;
            }

            // элементы на верхней границе
            j = mesh.Z.Count - 2;
            for (i = 1; i < mesh.X.Count - 2; i++)
            {
                current_index = i * (mesh.Z.Count - 1) + j;
                gamma_ij = gamma[i][j];

                neighbor_index = (i - 1) * (mesh.Z.Count - 1) + j - 1;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;

                neighbor_index = (i - 1) * (mesh.Z.Count - 1) + j;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;

                neighbor_index = i * (mesh.Z.Count - 1) + j - 1;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;

                neighbor_index = (i + 1) * (mesh.Z.Count - 1) + j - 1;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;

                neighbor_index = (i + 1) * (mesh.Z.Count - 1) + j;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;
            }

            // элементы на левой границе
            i = 0;
            for (j = 1; j < mesh.Z.Count - 2; j++)
            {
                current_index = i * (mesh.Z.Count - 1) + j;
                gamma_ij = gamma[i][j];

                neighbor_index = i * (mesh.Z.Count - 1) + j - 1;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;

                neighbor_index = i * (mesh.Z.Count - 1) + j + 1;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;

                neighbor_index = (i + 1) * (mesh.Z.Count - 1) + j - 1;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;

                neighbor_index = (i + 1) * (mesh.Z.Count - 1) + j;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;

                neighbor_index = (i + 1) * (mesh.Z.Count - 1) + j + 1;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;
            }

            // элементы на правой границе
            i = mesh.X.Count - 2;
            for (j = 1; j < mesh.Z.Count - 2; j++)
            {
                current_index = i * (mesh.Z.Count - 1) + j;
                gamma_ij = gamma[i][j];

                neighbor_index = (i - 1) * (mesh.Z.Count - 1) + j - 1;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;

                neighbor_index = (i - 1) * (mesh.Z.Count - 1) + j;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;

                neighbor_index = (i - 1) * (mesh.Z.Count - 1) + j + 1;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;


                neighbor_index = i * (mesh.Z.Count - 1) + j - 1;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;

                neighbor_index = i * (mesh.Z.Count - 1) + j + 1;
                difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
                functional += gamma_ij * difference * difference;
            }

            // левый нижний элемент   i = 0, j = 0
            i = 0;
            j = 0;
            current_index = i * (mesh.Z.Count - 1) + j;
            gamma_ij = gamma[i][j];

            neighbor_index = i * (mesh.Z.Count - 1) + j + 1;
            difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
            functional += gamma_ij * difference * difference;

            neighbor_index = (i + 1) * (mesh.Z.Count - 1) + j;
            difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
            functional += gamma_ij * difference * difference;

            neighbor_index = (i + 1) * (mesh.Z.Count - 1) + j + 1;
            difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
            functional += gamma_ij * difference * difference;


            // левый верхний элемент  i = 0, j = mesh.Z.Count - 2
            i = 0;
            j = mesh.Z.Count - 2;
            current_index = i * (mesh.Z.Count - 1) + j;
            gamma_ij = gamma[i][j];

            neighbor_index = i * (mesh.Z.Count - 1) + j - 1;
            difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
            functional += gamma_ij * difference * difference;

            neighbor_index = (i + 1) * (mesh.Z.Count - 1) + j - 1;
            difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
            functional += gamma_ij * difference * difference;

            neighbor_index = (i + 1) * (mesh.Z.Count - 1) + j;
            difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
            functional += gamma_ij * difference * difference;


            // правый нижний элемент  i = mesh.X.Count - 2, j = 0
            i = mesh.X.Count - 2;
            j = 0;
            current_index = i * (mesh.Z.Count - 1) + j;
            gamma_ij = gamma[i][j];

            neighbor_index = (i - 1) * (mesh.Z.Count - 1) + j;
            difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
            functional += gamma_ij * difference * difference;

            neighbor_index = (i - 1) * (mesh.Z.Count - 1) + j + 1;
            difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
            functional += gamma_ij * difference * difference;

            neighbor_index = i * (mesh.Z.Count - 1) + j + 1;
            difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
            functional += gamma_ij * difference * difference;


            // правый верхний элемент i = mesh.X.Count - 2, j = mesh.Z.Count - 2
            i = mesh.X.Count - 2;
            j = mesh.Z.Count - 2;
            current_index = i * (mesh.Z.Count - 1) + j;
            gamma_ij = gamma[i][j];

            neighbor_index = (i - 1) * (mesh.Z.Count - 1) + j - 1;
            difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
            functional += gamma_ij * difference * difference;

            neighbor_index = (i - 1) * (mesh.Z.Count - 1) + j;
            difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
            functional += gamma_ij * difference * difference;

            neighbor_index = i * (mesh.Z.Count - 1) + j - 1;
            difference = field_p_in_mesh_elements[current_index] - field_p_in_mesh_elements[neighbor_index];
            functional += gamma_ij * difference * difference;


            return functional;
        }


        public static void UpdateGamma(
            double[][] gamma,
            int rows, int cols,
            double[] field_p_in_mesh_elements,
            double gamma_mult,
            double allow_max_difference)
        {
            int current_index;
            double p_current, p_neighbor;
            double difference;
            double max_difference;

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    current_index = i * cols + j;
                    p_current = field_p_in_mesh_elements[current_index];
                    max_difference = 0.0;

                    if (i != 0)
                    {
                        if (j != 0)
                        {
                            p_neighbor = field_p_in_mesh_elements[(i - 1) * cols + j - 1];

                            if (Math.Abs(p_current) >= Math.Abs(p_neighbor))
                            {
                                difference = Math.Abs(p_current / p_neighbor);
                            }
                            else
                            {
                                difference = Math.Abs(p_neighbor / p_current);
                            }

                            max_difference = max_difference > difference ? max_difference : difference;
                        }


                        p_neighbor = field_p_in_mesh_elements[(i - 1) * cols + j];

                        if (Math.Abs(p_current) >= Math.Abs(p_neighbor))
                        {
                            difference = Math.Abs(p_current / p_neighbor);
                        }
                        else
                        {
                            difference = Math.Abs(p_neighbor / p_current);
                        }

                        max_difference = max_difference > difference ? max_difference : difference;

                        if (j != cols - 1)
                        {
                            p_neighbor = field_p_in_mesh_elements[(i - 1) * cols + j + 1];

                            if (Math.Abs(p_current) >= Math.Abs(p_neighbor))
                            {
                                difference = Math.Abs(p_current / p_neighbor);
                            }
                            else
                            {
                                difference = Math.Abs(p_neighbor / p_current);
                            }

                            max_difference = max_difference > difference ? max_difference : difference;
                        }
                    }

                    if (j != 0)
                    {
                        p_neighbor = field_p_in_mesh_elements[i * cols + j - 1];

                        if (Math.Abs(p_current) >= Math.Abs(p_neighbor))
                        {
                            difference = Math.Abs(p_current / p_neighbor);
                        }
                        else
                        {
                            difference = Math.Abs(p_neighbor / p_current);
                        }

                        max_difference = max_difference > difference ? max_difference : difference;
                    }

                    if (j != cols - 1)
                    {
                        p_neighbor = field_p_in_mesh_elements[i * cols + j + 1];

                        if (Math.Abs(p_current) >= Math.Abs(p_neighbor))
                        {
                            difference = Math.Abs(p_current / p_neighbor);
                        }
                        else
                        {
                            difference = Math.Abs(p_neighbor / p_current);
                        }

                        max_difference = max_difference > difference ? max_difference : difference;
                    }

                    if (i != rows - 1)
                    {
                        if (j != 0)
                        {
                            p_neighbor = field_p_in_mesh_elements[(i + 1) * cols + j - 1];

                            if (Math.Abs(p_current) >= Math.Abs(p_neighbor))
                            {
                                difference = Math.Abs(p_current / p_neighbor);
                            }
                            else
                            {
                                difference = Math.Abs(p_neighbor / p_current);
                            }

                            max_difference = max_difference > difference ? max_difference : difference;
                        }

                        p_neighbor = field_p_in_mesh_elements[(i + 1) * cols + j];

                        if (Math.Abs(p_current) >= Math.Abs(p_neighbor))
                        {
                            difference = Math.Abs(p_current / p_neighbor);
                        }
                        else
                        {
                            difference = Math.Abs(p_neighbor / p_current);
                        }

                        max_difference = max_difference > difference ? max_difference : difference;

                        if (j != cols - 1)
                        {
                            p_neighbor = field_p_in_mesh_elements[(i + 1) * cols + j + 1];

                            if (Math.Abs(p_current) >= Math.Abs(p_neighbor))
                            {
                                difference = Math.Abs(p_current / p_neighbor);
                            }
                            else
                            {
                                difference = Math.Abs(p_neighbor / p_current);
                            }

                            max_difference = max_difference > difference ? max_difference : difference;
                        }
                    }

                    if (max_difference > allow_max_difference)
                    {

                        gamma[i][j] *= gamma_mult;
                    }
                }
            }
        }
    }
}
