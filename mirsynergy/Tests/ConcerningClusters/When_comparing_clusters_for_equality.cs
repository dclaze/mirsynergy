using System.Collections.Generic;
using Xunit;
using Xunit.Should;

namespace mirsynergy.Tests.ConcerningClusters
{
    public class When_comparing_clusters_for_equality
    {
        [Fact]
        public void should_return_equal_if_contained_rna_are_the_same()
        {
            var cluster1 = new Cluster()
            {
                MicroRnaIndexes = new List<int>() { 4, 3, 2 }
            };
            var cluster2 = new Cluster()
            {
                MicroRnaIndexes = new List<int>() { 2, 3, 4 }
            };

            cluster1.Equals(cluster2).ShouldBeTrue();
        }

        [Fact]
        public void should_return_not_equal_if_contain_rna_are_different()
        {
            var cluster1 = new Cluster()
            {
                MicroRnaIndexes = new List<int>() { 1, 2, 3 }
            };
            var cluster2 = new Cluster()
            {
                MicroRnaIndexes = new List<int>() { 4, 5, 6 }
            };

            cluster1.Equals(cluster2).ShouldBeFalse();
        }

        [Fact]
        public void should_return_equal_if_both_clusters_have_no_indexes()
        {
            var cluster1 = new Cluster();
            var cluster2 = new Cluster();

            cluster1.Equals(cluster2).ShouldBeTrue();
        }

        [Fact]
        public void should_return_not_equal_if_one_cluster_has_indexes_and_another_cluster_has_none()
        {
            var cluster1 = new Cluster();
            var cluster2 = new Cluster() { MicroRnaIndexes = new List<int>() { 1 } };

            cluster1.Equals(cluster2).ShouldBeFalse();
        }
    }
}