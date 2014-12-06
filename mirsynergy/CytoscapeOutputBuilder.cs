using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace mirsynergy
{
    public class CytoscapeOutputBuilder
    {
        public CytoscapeOutputBuilder()
        {
            Random = new Random();
        }

        public Random Random { get; set; }

        public CytoscapeOutput Generate(List<Cluster> finalClusterAssignmentsFromLoadedScores, MatrixParser.MatrixResult loadedMicroRnaSynergyScores, MatrixParser.MatrixResult loadedGeneGeneSynergyScores)
        {
            var combinedScores = loadedMicroRnaSynergyScores.Matrix.DiagonalStack(loadedGeneGeneSynergyScores.Matrix);

            var microRnaNodes = loadedMicroRnaSynergyScores.ColumnLabels.Select((s, i) => new Elements.Node
            {
                Data = CreateNode(i, s, combinedScores)
            }).ToList();
            var nextStartingIndex = microRnaNodes.Count;
            var mRnaNodes = loadedGeneGeneSynergyScores.ColumnLabels.Select((s, i) => new Elements.Node()
            {
                Data = CreateNode((nextStartingIndex + i), s, combinedScores)
            });
            var nodes = microRnaNodes.Concat(mRnaNodes).ToList();

            var edges = new List<Elements.Edge>();
            foreach (Cluster cluster in finalClusterAssignmentsFromLoadedScores)
            {
                var clusterColor = GetRandomColor();
                var edgePairs = Permutations.GetAllPairs(cluster.MicroRnaIndexes);
                edges.AddRange(edgePairs.Select(tuple => new Elements.Edge()
                {
                    Data = new Elements.Edge.EdgeData()
                    {
                        Source = tuple.Item1.ToString(),
                        Target = tuple.Item2.ToString(),
                    },
                    Networks = new List<string>() { cluster.Id.ToString() }
                }));

                foreach (var node in cluster.MicroRnaIndexes)
                {
                    var matchinNode = nodes.Single(n => n.Data.Id == node.ToString());
                    matchinNode.Data.Score = SynergyCalculator.GetSynergyScore(combinedScores, cluster.MicroRnaIndexes);
                    matchinNode.Group = cluster.Id.ToString();
                    matchinNode.Data.Color = "#" + clusterColor.R.ToString("X2") + clusterColor.G.ToString("X2") + clusterColor.B.ToString("X2"); ;
                }
            }


            return new CytoscapeOutput()
            {
                Elements = new Elements()
                {
                    Nodes = nodes,
                    Edges = edges
                }
            };
        }

        private Color GetRandomColor()
        {
            var names = (KnownColor[])Enum.GetValues(typeof(KnownColor));
            var randomColorName = names[Random.Next(names.Length)];
            return Color.FromKnownColor(randomColorName);
        }

        private Elements.Node.NodeData CreateNode(int i, string s, Matrix<double> combinedScores)
        {
            return new Elements.Node.NodeData()
            {
                Id = i.ToString(),
                Name = s
            };
        }
    }
}