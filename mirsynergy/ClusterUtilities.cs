using System;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace mirsynergy
{
    public static class ClusterUtilities
    {
        public static double GetWeightsOfInternalEdges(Matrix<double> microRnaMicroRnaSynergyScores, List<int> microRnaIndexes)
        {
            var internalEdgePairs = GetInternalEdgePairs(microRnaIndexes);
            var totalWeightsOfInternalEdges = GetTotalWeight(microRnaMicroRnaSynergyScores, internalEdgePairs);
            return totalWeightsOfInternalEdges;
        }

        public static double GetTotalWeightsOfBoundaryEdges(Matrix<double> microRnaMicroRnaSynergyScores, List<int> microRnaIndexes)
        {
            var externalEdgePairs = GetExternalEdgePairs(microRnaMicroRnaSynergyScores, microRnaIndexes);
            var totalWeightsOfBoundaryEdges = GetTotalWeight(microRnaMicroRnaSynergyScores, externalEdgePairs);
            return totalWeightsOfBoundaryEdges;
        }

        public static double OverlapScore(Cluster leftCluster, Cluster rightCluster)
        {
            var intersectionCardinality = Math.Pow(leftCluster.MicroRnaIndexes.Intersect(rightCluster.MicroRnaIndexes).Count(), 2);
            var unionCardinality = leftCluster.MicroRnaIndexes.Union(rightCluster.MicroRnaIndexes).Count();
            return intersectionCardinality / unionCardinality;
        }

        public static double Density1(Cluster cluster, Matrix<double> microRnaMicroRnaSynergyScores)
        {
            return (2 * GetWeightsOfInternalEdges(microRnaMicroRnaSynergyScores, cluster.MicroRnaIndexes)) / (cluster.MicroRnaIndexes.Count * (cluster.MicroRnaIndexes.Count - 1));
        }

        public static double Density2(Cluster cluster, Matrix<double> microRnaMicroRnaSynergyScores, Matrix<double> mRnamRnaSynergyScores)
        {
            return (2 * GetWeightsOfInternalEdges(microRnaMicroRnaSynergyScores, cluster.MicroRnaIndexes)) / (cluster.MicroRnaIndexes.Count * (cluster.MicroRnaIndexes.Count - 1));
        }

        private static double GetTotalWeight(Matrix<double> microRnaMicroRnaSynergyScores, IEnumerable<Tuple<int, int>> allPairs)
        {
            return allPairs.Sum(tuple => microRnaMicroRnaSynergyScores[tuple.Item1, tuple.Item2]);
        }

        private static IEnumerable<Tuple<int, int>> GetInternalEdgePairs(List<int> microRnaIndexes)
        {
            return Permutations.GetAllPairs(microRnaIndexes);
        }
        private static IEnumerable<Tuple<int, int>> GetExternalEdgePairs(Matrix<double> microRnaMicroRnaSynergyScores, ICollection<int> microRnaIndexes)
        {
            var externalEdgePairs = new List<Tuple<int, int>>();
            for (var externalIndex = 0; externalIndex < microRnaMicroRnaSynergyScores.RowCount; externalIndex++)
            {
                if (microRnaIndexes.Contains(externalIndex))
                    continue;
                externalEdgePairs.AddRange(
                    microRnaIndexes.Select(internalIndexes => new Tuple<int, int>(internalIndexes, externalIndex)));
            }

            return externalEdgePairs;
        }
    }
}