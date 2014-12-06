using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;

namespace mirsynergy
{
    public static class SynergyCalculator
    {
        public static double GetSynergyScore(Matrix<double> microRnaMicroRnaSynergyScores, List<int> microRnaIndexes)
        {
            var totalWeightsOfInternalEdges = ClusterUtilities.GetWeightsOfInternalEdges(microRnaMicroRnaSynergyScores, microRnaIndexes);

            var totalWeightsOfBoundaryEdges = ClusterUtilities.GetTotalWeightsOfBoundaryEdges(microRnaMicroRnaSynergyScores, microRnaIndexes);

            var penaltyScoreForFormingCluster = GetPenaltyScore(microRnaIndexes);


            var numerator = totalWeightsOfInternalEdges;
            var denomonator = (totalWeightsOfInternalEdges + totalWeightsOfBoundaryEdges + penaltyScoreForFormingCluster);
            if (numerator == 0 || denomonator == 0)
                return 0;

            return numerator / denomonator;
        }

        public static double GetPenaltyScore(IReadOnlyCollection<int> microRnaIndexes)
        {
            return 2 * microRnaIndexes.Count;
        }
    }
}