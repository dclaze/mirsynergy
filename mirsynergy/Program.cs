using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;


namespace mirsynergy
{
    class Program
    {
        static void Main(string[] args)
        {
            var expressions = new List<string>();
            var microRnas = new List<string>();

            var microRnaExpressionsMatrix = DenseMatrix.OfArray(new double[,]
            {
                {1.3, 1.1, 1.1},
                {1.2, 2.1, 7.0},
                {4.2, 3.2, 2.0}
            });
            var lassoScoringMatrix = CalculateLassoAndGetMmiw();
            var microRnaMicroRnaSynergyScoresMatrix = CalculateMiRnaMiRnaSynergisticScores(microRnaExpressionsMatrix, lassoScoringMatrix);

            //            var geneGeneInteraction = CalculateGeneGeneInteraction();
            var stage1Clusters = GetStage1Clusters(microRnaMicroRnaSynergyScoresMatrix, microRnas);
            var finalClusterAssignments = GetFinalClusterAssignments(null);
        }

        private static object GetFinalClusterAssignments(Matrix microRnaMicroRnaSynergyScores)
        {
            var clusters = new List<Cluster>();
            var availableMicroRnas = microRnaMicroRnaSynergyScores.EnumerateRowsIndexed().Select(tuple => tuple.Item1).ToList();

            while (availableMicroRnas.Any())
            {
                var microRnaIndexWithMaxSynergyScore = GetMicroRnaWithMaxTotalSynergyScores(microRnaMicroRnaSynergyScores, availableMicroRnas);
                var currentCluster = new Cluster(microRnaIndexWithMaxSynergyScore);
                var previousCluster = new Cluster();
                while (!previousCluster.Equals(currentCluster))
                {
                    previousCluster.MicroRnaIndexes = currentCluster.MicroRnaIndexes.ToList();

                    var bestmiRnaNeighbor = ChooseBestNeighboringMicroRna(microRnaMicroRnaSynergyScores, previousCluster);
                    var bestmiRnaToRemove = ChooseBestMicroRnaToRemove(microRnaMicroRnaSynergyScores, previousCluster);

                    var synergyOfPreviousClusterAndBestNeighbor = GetSynergyScore(microRnaMicroRnaSynergyScores,
                        previousCluster.MicroRnaIndexes.Concat(new[] { bestmiRnaNeighbor }).ToList());
                    var synergyOfPreviousCluster = GetSynergyScore(microRnaMicroRnaSynergyScores,
                        previousCluster.MicroRnaIndexes);
                    var synergyOfPreviousClusterWitBestRnaRemoved = GetSynergyScore(microRnaMicroRnaSynergyScores,
                        previousCluster.MicroRnaIndexes.Except(new[] { bestmiRnaToRemove }).ToList());

                    if (synergyOfPreviousClusterAndBestNeighbor >
                        Math.Max(synergyOfPreviousClusterWitBestRnaRemoved, synergyOfPreviousCluster))
                    {
                        currentCluster.ReplaceWith(previousCluster.MicroRnaIndexes.Union(new[] { bestmiRnaNeighbor }));
                    }
                    else if (synergyOfPreviousClusterWitBestRnaRemoved >
                             Math.Max(synergyOfPreviousClusterAndBestNeighbor, synergyOfPreviousCluster))
                    {
                        currentCluster.ReplaceWith(previousCluster.MicroRnaIndexes.Except(new[] { bestmiRnaToRemove }));
                    }
                }
                clusters = clusters.Union(new[] { previousCluster }).ToList();
                availableMicroRnas = availableMicroRnas.Except(new[] { microRnaIndexWithMaxSynergyScore }).ToList();
            }

            return clusters;
        }

        private static List<Cluster> GetStage1Clusters(Matrix microRnaMicroRnaSynergyScores, IEnumerable<string> microRnas)
        {
            var mirmClusters = GrowMirmsByOverlappingNeighborhoodExpansion(microRnaMicroRnaSynergyScores);
            return MergeMirmsByBreadthFirstSearch(mirmClusters);
        }

        private static List<Cluster> MergeMirmsByBreadthFirstSearch(List<Cluster> mirmClusters)
        {
            var mergedMirms = mirmClusters.ToList();
            var complementMergedMirms = new List<Cluster>();
            var threshold = 0.8;

            foreach (var currentCluster in mirmClusters)
            {
                var closelyConnectedClusters = new List<Cluster> { currentCluster };
                mergedMirms = mergedMirms.Except(new[] { currentCluster }).ToList();

                foreach (var otherCluster in mirmClusters.Except(new[] { currentCluster }))
                {
                    if (ClusterUtilities.OverlapScore(currentCluster, otherCluster) >= threshold)
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

        private static List<Cluster> GrowMirmsByOverlappingNeighborhoodExpansion(Matrix microRnaMicroRnaSynergyScores)
        {
            var clusters = new List<Cluster>();
            var availableMicroRnas = microRnaMicroRnaSynergyScores.EnumerateRowsIndexed().Select(tuple => tuple.Item1).ToList();

            while (availableMicroRnas.Any())
            {
                var microRnaIndexWithMaxSynergyScore = GetMicroRnaWithMaxTotalSynergyScores(microRnaMicroRnaSynergyScores, availableMicroRnas);
                var previousCluster = GetNextClusterUsingNeighboringExpansion(microRnaMicroRnaSynergyScores, microRnaIndexWithMaxSynergyScore);
                clusters = clusters.Union(new[] { previousCluster }).ToList();
                availableMicroRnas = availableMicroRnas.Except(new[] { microRnaIndexWithMaxSynergyScore }).ToList();
            }

            return clusters;
        }

        private static Cluster GetNextClusterUsingNeighboringExpansion(Matrix microRnaMicroRnaSynergyScores, int nextMicroRnaIndexSeed)
        {
            var currentCluster = new Cluster(nextMicroRnaIndexSeed);
            var previousCluster = new Cluster();
            while (!previousCluster.Equals(currentCluster))
            {
                previousCluster.MicroRnaIndexes = currentCluster.MicroRnaIndexes.ToList();

                var bestmiRnaNeighbor = ChooseBestNeighboringMicroRna(microRnaMicroRnaSynergyScores, previousCluster);
                var bestmiRnaToRemove = ChooseBestMicroRnaToRemove(microRnaMicroRnaSynergyScores, previousCluster);

                var synergyOfPreviousClusterAndBestNeighbor = GetSynergyScore(microRnaMicroRnaSynergyScores,
                    previousCluster.MicroRnaIndexes.Concat(new[] { bestmiRnaNeighbor }).ToList());
                var synergyOfPreviousCluster = GetSynergyScore(microRnaMicroRnaSynergyScores, previousCluster.MicroRnaIndexes);
                var synergyOfPreviousClusterWitBestRnaRemoved = GetSynergyScore(microRnaMicroRnaSynergyScores,
                    previousCluster.MicroRnaIndexes.Except(new[] { bestmiRnaToRemove }).ToList());

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

        private static int ChooseBestMicroRnaToRemove(Matrix microRnaMicroRnaSynergyScores, Cluster previousCluster)
        {
            var clustersWithOneRemoved = previousCluster.MicroRnaIndexes.ToDictionary(microRnaIndex => microRnaIndex,
                microRnaIndex => previousCluster.MicroRnaIndexes.Except(new[] { microRnaIndex }).ToList());

            return clustersWithOneRemoved.OrderByDescending(pair => GetSynergyScore(microRnaMicroRnaSynergyScores, pair.Value)).First().Key;
        }

        private static int ChooseBestNeighboringMicroRna(Matrix microRnaMicroRnaSynergyScores, Cluster previousCluster)
        {
            var indexesNotIncludedInCluster = GetOtherIndexes(microRnaMicroRnaSynergyScores, previousCluster.MicroRnaIndexes);
            //            if(indexesNotIncludedInCluster.Any())
            return indexesNotIncludedInCluster.OrderByDescending(i => GetSynergyScore(microRnaMicroRnaSynergyScores, previousCluster.MicroRnaIndexes.Concat(new[] { i }).ToList())).First();
        }

        private static double GetSynergyScore(Matrix microRnaMicroRnaSynergyScores, List<int> microRnaIndexes)
        {
            //TODO Filter out zero weights
            var totalWeightsOfInternalEdges = ClusterUtilities.GetWeightsOfInternalEdges(microRnaMicroRnaSynergyScores, microRnaIndexes);

            var totalWeightsOfBoundaryEdges = ClusterUtilities.GetTotalWeightsOfBoundaryEdges(microRnaMicroRnaSynergyScores, microRnaIndexes);

            var penaltyScoreForFormingCluster = GetPenaltyScore(microRnaIndexes);


            var numerator = totalWeightsOfInternalEdges;
            var denomonator = (totalWeightsOfInternalEdges + totalWeightsOfBoundaryEdges + penaltyScoreForFormingCluster);
            if (numerator == 0 || denomonator == 0)
                return 0;

            return numerator / denomonator;
        }




        private static double GetPenaltyScore(IReadOnlyCollection<int> microRnaIndexes)
        {
            return 2 * microRnaIndexes.Count;
        }





        private static List<int> GetOtherIndexes(Matrix microRnaMicroRnaSynergyScores, List<int> microRnaIndexes)
        {
            var excludedIndexes = new List<int>();

            for (var externalIndex = 0; externalIndex < microRnaMicroRnaSynergyScores.RowCount; externalIndex++)
            {
                if (microRnaIndexes.Contains(externalIndex))
                    continue;
                excludedIndexes.Add(externalIndex);
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

        private static Matrix<double> CalculateGeneGeneInteraction()
        {
            throw new NotImplementedException();
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

        //miRNA-mRNA interaction weight MMIW
        private static Matrix<double> CalculateLassoAndGetMmiw()
        {
            return Matrix<double>.Build.Random(500, 500);
        }
    }

    public class MicroRNA
    {
        public string Id { get; set; }
    }

    public class SquareMatrixOperationAttemptedOnNonSquareMatrix : Exception
    {
    }
}
