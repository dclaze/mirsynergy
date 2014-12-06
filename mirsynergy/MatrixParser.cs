using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace mirsynergy
{
    public static class MatrixParser
    {
        public static MatrixResult ParseFromFile(string fileName)
        {
            var values = File.ReadAllText(fileName).Split(new string[] { "\n", "\r\n" }, StringSplitOptions.RemoveEmptyEntries).Select(s => s.Split(',')).ToList();
            var columnLabels = values.First().Skip(1).ToList();
            values = values.Skip(1).ToList();
            var rowLabels = values.Select(strings => strings.First()).ToList();
            values = values.Select(strings => strings.Skip(1).ToArray()).ToList();

            var rowCount = rowLabels.Count();
            var columnCount = columnLabels.Count();

            var matrix = new double[rowCount, columnCount];
            for (var row = 0; row < rowCount; row++)
            {
                for (var column = 0; column < columnCount; column++)
                {
                    var s = values[row][column];
                    matrix[row, column] = Double.Parse(s);
                }
            }
            return new MatrixResult()
            {
                Matrix = DenseMatrix.OfArray(matrix),
                RowLabels = rowLabels,
                ColumnLabels = columnLabels
            };
        }

        public class MatrixResult
        {
            public Matrix<double> Matrix { get; set; }
            public List<string> RowLabels { get; set; }
            public List<string> ColumnLabels { get; set; }
        }
    }
}