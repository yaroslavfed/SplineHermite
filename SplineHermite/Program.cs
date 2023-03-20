using System.Diagnostics;
using System.Drawing;
using static System.Math;

namespace SplineHermite
{
    internal class Program
    {
        static List<double> points = new List<double>();    // известные точки
        static List<double> x = new List<double>();         // весь х
        static List<double> Hfx = new List<double>();       // значение полинома в x

        static double globalStep = 0.001;                   // шаг для вывода графика

        static private void ReadFile(string path)
        {
            using (StreamReader sr = new StreamReader(path))
            {
                foreach (var number in sr.ReadLine()!.Split(' '))
                {
                    points.Add(double.Parse(number));
                }
            };
        }

        #region DrawPlot
        static private void DrawPlot()
        {
            var result = x.Concat(Hfx);

            using Process myProcess = new Process();
            myProcess.StartInfo.FileName = "python";
            myProcess.StartInfo.Arguments = @"script.py";
            myProcess.StartInfo.UseShellExecute = false;
            myProcess.StartInfo.RedirectStandardInput = true;
            myProcess.StartInfo.RedirectStandardOutput = false;
            myProcess.Start();

            string outputPath = Path.Combine(Directory.GetCurrentDirectory(), "test.txt");
            using (StreamWriter sw = new StreamWriter(outputPath))
            {
                foreach (var point in points)
                    sw.Write(point + " ");
                sw.WriteLine();

                foreach (var point in points)
                    sw.Write(function(point) + " ");
                sw.WriteLine();

                foreach (var point in x)
                    sw.Write(point + " ");
                sw.WriteLine();

                foreach (var point in Hfx)
                    sw.Write(point + " ");
                sw.WriteLine();
            };
        }
        #endregion DrawPlot

        static private double function(double x) => 1 / (1 + 25 * Pow(x, 2));       // Pow(E, Sin(PI * x))   Pow(x, 3)   Cos(Sin(Pow(E, x)))

        #region PiecewiseCubicHermite
        static private double eps(double x, double xi, double h) => (x - xi) / h;
        static private double dx(double fiplus, double fi, double fimin, double hi, double himin) => (1/2) * (((fiplus - fi)/hi) + ((fi - fimin)/himin));
        static private double dx0(double fiplusplus, double fiplus, double fi, double hi, double hiplus) => (1/2) * (-((3*hi + 2*hiplus)/(hi*(hi + hiplus)))*fi + ((hi + 2*hiplus)/(hi*hiplus))*fiplus - (hi/((hi + hiplus)*hiplus))*fiplusplus);
        static private double dxn(double fiminmin, double fimin, double fi, double himinmin, double himin) => (1/2) * (((himin)/(himinmin*(himinmin + himin)))*fiminmin - ((2*himinmin + himin)/(himin*himinmin))*fimin + ((3*himin + 2*himinmin)/(himin*(himinmin + himin)))*fi);

        static private double BasicFunctions(int i, double x, double xi, double xiplus, double h)
        {
            switch (i)
            {
                case 0:
                    return 1 - (3 * Pow(eps(x, xi, h), 2)) + (2 * Pow(eps(x, xi, h), 3));
                case 1:
                    return h * (eps(x, xi, h) - 2 * Pow(eps(x, xi, h), 2) + Pow(eps(x, xi, h), 3));
                case 2:
                    return (3 * Pow(eps(x, xi, h), 2)) - (2 * Pow(eps(x, xi, h), 3));
                case 3:
                    return h * ( - Pow(eps(x, xi, h), 2) + Pow(eps(x, xi, h), 3));
                default:
                    return 0;
            }
        }

        static private void PiecewiseCubicHermite()
        {
            Console.Write("Узлы: ");
            foreach (double point in points)
                Console.Write("{0} ", point);
            Console.WriteLine();

            int n = points.Count;

            double result = 0;
            double derivative = 1;
            double derivativePlus = 1;

            for (int p = 0; p < n - 1; p++)
            {
                double step = points[p + 1] - points[p];

                for (double k = points[p]; k <= points[p + 1]; k += globalStep)
                {
                    x.Add(k);
                    if(p == 0)
                    {
                        derivative = dx0(function(points[p + 2]), function(points[p + 1]), function(points[p]), step, points[p + 2] - points[p + 1]);
                    }
                    else
                    {
                        derivative = dx(function(points[p + 1]), function(points[p]), function(points[p - 1]), step, points[p] - points[p - 1]);
                    }
                    if (p == n - 2)
                    {
                        derivativePlus = dxn(function(points[n-1 - 2]), function(points[n-1 - 1]), function(points[n-1]), points[n-1 - 1] - points[n-1 - 2], points[n-1] - points[n-1 - 1]);
                    }
                    else
                    {
                        derivativePlus = dx(function(points[p + 1 + 1]), function(points[p + 1]), function(points[p + 1 - 1]), points[p + 1 + 1] - points[p + 1], points[p + 1] - points[p + 1 - 1]);
                    }
                    

                    result = function(points[p]) * BasicFunctions(0, k, points[p], points[p + 1], step) +
                        derivative * BasicFunctions(1, k, points[p], points[p + 1], step) +
                        function(points[p + 1]) * BasicFunctions(2, k, points[p], points[p + 1], step) +
                        derivativePlus * BasicFunctions(3, k, points[p], points[p + 1], step);

                    Console.WriteLine("x: {0}\ty: {1}", k, result);
                    Hfx.Add(result);
                }
            }
            result = function(points[n - 1]);
            Console.WriteLine("x: {0}\ty: {1}", points[n - 1], result);
            x.Add(points[n - 1]);
            Hfx.Add(result);
        }
        #endregion PiecewiseCubicHermite

        static void Main(string[] args)
        {
            string inputPath = Path.Combine(Directory.GetCurrentDirectory(), "points.txt");
            ReadFile(inputPath);
            PiecewiseCubicHermite();
            DrawPlot();
        }
    }
}