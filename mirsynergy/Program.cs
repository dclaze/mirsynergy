using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;


namespace mirsynergy
{
    class Program
    {
        private const double OverlapThreshold = 0.8;
        private const double Density2Threshold = 5e-3;
        private const double Density1Threshold = 1e-2;

        static void Main(string[] args)
        {
            var microRnaSynergyScores = MatrixParser.ParseFromFile("toy_modules_W.csv");
            var geneGeneSynergyScores = MatrixParser.ParseFromFile("toy_modules_H.csv");
            var stage1ClustersFromLoadScores = GetStage1Clusters(microRnaSynergyScores.Matrix);
            var finalClusterAssignmentsFromLoadedScores = GetFinalClusterAssignments(stage1ClustersFromLoadScores, microRnaSynergyScores.Matrix, geneGeneSynergyScores.Matrix);
            var cytoscapeBuilder = new CytoscapeOutputBuilder();
            var ctyoscapeOutput = cytoscapeBuilder.Generate(finalClusterAssignmentsFromLoadedScores, microRnaSynergyScores, geneGeneSynergyScores);
            JsonFileWriter.WriteToFile(ctyoscapeOutput, "output.json");
        }

        private static List<Cluster> GetStage1Clusters(Matrix<double> microRnaMicroRnaSynergyScores)
        {
            var mirmClusters = GrowMirmsByOverlappingNeighborhoodExpansion(microRnaMicroRnaSynergyScores);
            var mergedMirmClusters = MergeMirmsByBreadthFirstSearch(mirmClusters);
            return mergedMirmClusters.Where(cluster => ClusterUtilities.Density1(cluster, microRnaMicroRnaSynergyScores) >= Density1Threshold).ToList();
        }

        private static IEnumerable<Cluster> MergeMirmsByBreadthFirstSearch(List<Cluster> mirmClusters)
        {
            var mergedMirms = mirmClusters.ToList();
            var complementMergedMirms = new List<Cluster>();
            foreach (var currentCluster in mirmClusters)
            {
                var closelyConnectedClusters = new List<Cluster> { currentCluster };
                mergedMirms = mergedMirms.Except(new[] { currentCluster }).ToList();

                foreach (var otherCluster in mirmClusters.Except(new[] { currentCluster }))
                {
                    if (ClusterUtilities.OverlapScore(currentCluster, otherCluster) >= OverlapThreshold)
                    {
                        closelyConnectedClusters = closelyConnectedClusters.Union(new[] { otherCluster }).ToList();
                        mergedMirms = mergedMirms.Except(new[] { otherCluster }).ToList();
                    }
                }

                var complementCloselyConnectedClusters = new List<Cluster>();
                foreach (var closelyConnectedCluster in closelyConnectedClusters)
                {
                    complementCloselyConnectedClusters =
                        new List<Cluster>(complementCloselyConnectedClusters.Union(new[] { closelyConnectedCluster })).ToList();
                }

                complementMergedMirms = complementMergedMirms.Union(complementCloselyConnectedClusters).ToList();
            }

            return complementMergedMirms;
        }

        private static List<Cluster> GetFinalClusterAssignments(IEnumerable<Cluster> stage1Clusters, Matrix<double> geneGeneInteractionSynergyScores, Matrix<double> microRnaSynergyScores)
        {
            var clusters = new List<Cluster>();
            var combinedMicroRnaAndmRnaSynergyScoresMatrix = microRnaSynergyScores.DiagonalStack(geneGeneInteractionSynergyScores);

            foreach (var cluster in stage1Clusters)
            {
                var previousCluster = GetNextClusterUsingNeighboringExpansion(combinedMicroRnaAndmRnaSynergyScoresMatrix, cluster, geneGeneInteractionSynergyScores.ColumnCount);
                clusters = clusters.Union(new[] { previousCluster }).ToList();
            }

            return clusters.Where(cluster => ClusterUtilities.Density2(cluster, microRnaSynergyScores, geneGeneInteractionSynergyScores)>=Density2Threshold).ToList();
        }

        private static List<Cluster> GrowMirmsByOverlappingNeighborhoodExpansion(Matrix<double> microRnaMicroRnaSynergyScores)
        {
            var clusters = new List<Cluster>();
            var availableMicroRnas = microRnaMicroRnaSynergyScores.EnumerateRowsIndexed().Select(tuple => tuple.Item1).ToList();

            while (availableMicroRnas.Any())
            {
                var microRnaIndexWithMaxSynergyScore = GetMicroRnaWithMaxTotalSynergyScores(microRnaMicroRnaSynergyScores, availableMicroRnas);
                var previousCluster = GetNextClusterUsingNeighboringExpansion(microRnaMicroRnaSynergyScores, new Cluster(microRnaIndexWithMaxSynergyScore));
                clusters = clusters.Union(new[] { previousCluster }).ToList();
                availableMicroRnas.Remove(microRnaIndexWithMaxSynergyScore);
            }

            return clusters;
        }

        private static Cluster GetNextClusterUsingNeighboringExpansion(Matrix<double> microRnaMicroRnaSynergyScores, Cluster currentCluster, int minumumIndex = 0)
        {
            var previousCluster = new Cluster();
            while (!previousCluster.Equals(currentCluster))
            {
                previousCluster.MicroRnaIndexes = currentCluster.MicroRnaIndexes.ToList();

                var bestmiRnaNeighbor = ChooseBestNeighboringMicroRna(microRnaMicroRnaSynergyScores, previousCluster, minumumIndex);
                var synergyOfPreviousClusterAndBestNeighbor = bestmiRnaNeighbor > 0 ? SynergyCalculator.GetSynergyScore(microRnaMicroRnaSynergyScores,
                    previousCluster.MicroRnaIndexes.Concat(new[] { bestmiRnaNeighbor }).ToList()) : 0;

                var synergyOfPreviousCluster = SynergyCalculator.GetSynergyScore(microRnaMicroRnaSynergyScores, previousCluster.MicroRnaIndexes);

                var bestmiRnaToRemove = ChooseBestIndexToRemove(microRnaMicroRnaSynergyScores, previousCluster, minumumIndex);
                var synergyOfPreviousClusterWitBestRnaRemoved = bestmiRnaToRemove > 0 ? SynergyCalculator.GetSynergyScore(microRnaMicroRnaSynergyScores,
                    previousCluster.MicroRnaIndexes.Except(new[] { bestmiRnaToRemove }).ToList()) : 0;

                if (synergyOfPreviousClusterAndBestNeighbor >
                    Math.Max(synergyOfPreviousClusterWitBestRnaRemoved, synergyOfPreviousCluster))
                {
                    currentCluster.ReplaceWith(previousCluster.MicroRnaIndexes.Union(new[] { bestmiRnaNeighbor }));
                }
                else if (synergyOfPreviousClusterWitBestRnaRemoved > Math.Max(synergyOfPreviousClusterAndBestNeighbor, synergyOfPreviousCluster))
                {
                    currentCluster.ReplaceWith(previousCluster.MicroRnaIndexes.Except(new[] { bestmiRnaToRemove }));
                }
            }
            return previousCluster;
        }

        private static int ChooseBestIndexToRemove(Matrix<double> microRnaMicroRnaSynergyScores, Cluster previousCluster, int minumumIndex)
        {
            var clustersWithOneRemoved = previousCluster.MicroRnaIndexes.Where(i => i >= minumumIndex).ToDictionary(microRnaIndex => microRnaIndex,
                microRnaIndex => previousCluster.MicroRnaIndexes.Except(new[] { microRnaIndex }).ToList());
            if (!clustersWithOneRemoved.Any())
                return -1;
            return clustersWithOneRemoved.OrderByDescending(pair => SynergyCalculator.GetSynergyScore(microRnaMicroRnaSynergyScores, pair.Value)).First().Key;
        }

        private static int ChooseBestNeighboringMicroRna(Matrix<double> microRnaMicroRnaSynergyScores, Cluster previousCluster, int minumumIndex)
        {
            var indexesNotIncludedInCluster = GetNeighboringIndexes(microRnaMicroRnaSynergyScores, previousCluster.MicroRnaIndexes, minumumIndex);
            if (!indexesNotIncludedInCluster.Any())
                return -1;
            return indexesNotIncludedInCluster.OrderByDescending(i => SynergyCalculator.GetSynergyScore(microRnaMicroRnaSynergyScores, previousCluster.MicroRnaIndexes.Concat(new[] { i }).ToList())).First();
        }



        private static List<int> GetNeighboringIndexes(Matrix<double> microRnaMicroRnaSynergyScores, List<int> microRnaIndexes, int minumumIndex)
        {
            var excludedIndexes = new List<int>();

            for (var externalIndex = minumumIndex; externalIndex < microRnaMicroRnaSynergyScores.RowCount; externalIndex++)
            {
                if (microRnaIndexes.Contains(externalIndex))
                    continue;
                excludedIndexes.AddRange(from microRnaIndex in microRnaIndexes
                                         where microRnaMicroRnaSynergyScores[microRnaIndex, externalIndex] != 0
                                         select externalIndex);
            }

            return excludedIndexes;
        }

        private static int GetMicroRnaWithMaxTotalSynergyScores(Matrix<double> microRnaMicroRnaSynergyScores, List<int> availableMicroRnas)
        {
            var indexOfMicroRnaWithMaxSynergyScore = 0;
            var currentMax = Double.MinValue;
            foreach (var column in microRnaMicroRnaSynergyScores.EnumerateColumnsIndexed())
            {
                if (!availableMicroRnas.Contains(column.Item1)) continue;

                var columnSum = column.Item2.Sum();

                if (!(columnSum > currentMax)) continue;

                indexOfMicroRnaWithMaxSynergyScore = column.Item1;
                currentMax = columnSum;
            }

            return indexOfMicroRnaWithMaxSynergyScore;
        }

        private static Matrix CalculateMiRnaMiRnaSynergisticScores(Matrix<double> microRnaExpressions, Matrix<double> lassoScoringMatrix)
        {
            var synergyisticScores = DenseMatrix.Create(microRnaExpressions.RowCount, microRnaExpressions.ColumnCount, 0);

            for (var rowIndex = 0; rowIndex < microRnaExpressions.RowCount; rowIndex++)
            {
                for (var columnIndex = 0; columnIndex < microRnaExpressions.ColumnCount; columnIndex++)
                {
                    synergyisticScores[rowIndex, columnIndex] = GetSynergyScore(lassoScoringMatrix, rowIndex, columnIndex);
                }
            }
            return synergyisticScores;
        }

        private static double GetSynergyScore(Matrix<double> lassoScoringMatrix, int rowIndex, int columnIndex)
        {
            if (!IsSquare(lassoScoringMatrix))
                throw new SquareMatrixOperationAttemptedOnNonSquareMatrix();
            double sumOfLassoScoresForRowAndColumn = 0;
            double sumOfLassoScoreForRow = 0;
            double sumOfLassoScoreForColumn = 0;

            for (var i = 0; i < lassoScoringMatrix.RowCount; i++)
            {
                sumOfLassoScoresForRowAndColumn += lassoScoringMatrix[i, rowIndex] * lassoScoringMatrix[i, columnIndex];
                sumOfLassoScoreForRow += lassoScoringMatrix[i, rowIndex];
                sumOfLassoScoreForColumn += lassoScoringMatrix[i, columnIndex];
            }

            return sumOfLassoScoresForRowAndColumn / Math.Min(sumOfLassoScoreForRow, sumOfLassoScoreForColumn);
        }

        private static bool IsSquare(Matrix<double> lassoScoringMatrix)
        {
            return lassoScoringMatrix.RowCount == lassoScoringMatrix.ColumnCount;
        }
    }

    public class SquareMatrixOperationAttemptedOnNonSquareMatrix : Exception
    {
    }
}
