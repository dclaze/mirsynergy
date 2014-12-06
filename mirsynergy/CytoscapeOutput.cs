using System.Collections.Generic;

namespace mirsynergy
{
    public class CytoscapeOutput
    {
        public Elements Elements { get; set; }
    }

    public class Elements
    {
        public List<Node> Nodes { get; set; }
        public List<Edge> Edges { get; set; }
        public class Node
        {
            public NodeData Data { get; set; }
            public class NodeData
            {
                public string Id { get; set; }
                public string Name { get; set; }
                public double Score { get; set; }
                public string Color = "white";
            }
            public string Group { get; set; }
            
        }

        public class Edge
        {
            public EdgeData Data { get; set; }
            public List<string> Networks { get; set; }
            public class EdgeData
            {
                public string Id { get; set; }
                public string Source { get; set; }
                public string Target { get; set; }
            }
        }
    }
}