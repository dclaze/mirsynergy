using System;
using System.Collections.Generic;
using System.Linq;

namespace mirsynergy
{
    public class Cluster
    {
        public List<int> MicroRnaIndexes { get; set; }

        public Cluster()
        {
            MicroRnaIndexes = new List<int>();
        }

        public Cluster(params int[] microRnaIndexWithMaxSynergyScore)
        {
            MicroRnaIndexes = new List<int>(microRnaIndexWithMaxSynergyScore);
        }

        public void ReplaceWith(IEnumerable<int> newMicroRnaIndexes)
        {
            MicroRnaIndexes = newMicroRnaIndexes.ToList();
        }

        protected bool Equals(Cluster other)
        {
            return Equals(MicroRnaIndexes, other.MicroRnaIndexes);
        }

        public override bool Equals(Object obj)
        {
            var otherCluster = (Cluster) obj;
            return MicroRnaIndexes.Count == otherCluster.MicroRnaIndexes.Count && !MicroRnaIndexes.Except(otherCluster.MicroRnaIndexes).Any();
        }

        public override int GetHashCode()
        {
            return MicroRnaIndexes.GetHashCode();
        }
    }
}