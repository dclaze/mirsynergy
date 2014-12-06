using System;
using Xunit;
using Xunit.Should;

namespace mirsynergy.Tests.ConcerningReadMatricies
{
    public class When_reading_a_matrix_from_a_file
    {
        [Fact]
        public void it_should_properly_serialize()
        {
            var matrix = MatrixParser.ParseFromFile("Z.csv");
            matrix.Matrix.ShouldNotBeNull();
        }
    }
}