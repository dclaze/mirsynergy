using System;
using Xunit;
using Xunit.Should;

namespace mirsynergy.Tests
{
    public class When_generating_all_pairs
    {
        [Fact]
        public void should_generate_correct_number_of_pairs()
        {
            var pairs = Permutations.GetAllPairs(new[] { "A", "B", "C", "D", "E" });
            pairs.Count.ShouldBe(10);
        }

        [Fact]
        public void should_generate_correct_pairs()
        {
            var pairs = Permutations.GetAllPairs(new[] { "A", "B", "C" });

            pairs.ShouldContain(new Tuple<string, string>("A", "B"));
            pairs.ShouldContain(new Tuple<string, string>("A", "C"));
            pairs.ShouldContain(new Tuple<string, string>("B", "C"));
            pairs.Count.ShouldBe(3);
        }

        [Fact]
        public void should_generate_one_pair_given_two_items()
        {
            var pairs = Permutations.GetAllPairs(new[] { "A", "B" });
            pairs.ShouldContain(new Tuple<string, string>("A", "B"));
            pairs.Count.ShouldBe(1);
        }

        [Fact]
        public void should_generate_no_pair_given_one_item()
        {
            var pairs = Permutations.GetAllPairs(new[] { "F" });
            pairs.ShouldBeEmpty();
        }

        [Fact]
        public void should_generate_no_pair_given_no_item()
        {
            var strings = new string[0];
            var pairs = Permutations.GetAllPairs(strings);
            pairs.ShouldBeEmpty();
        }
    }
}