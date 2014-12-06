using System;
using System.Collections.Generic;
using System.Linq;

namespace mirsynergy
{
    public static class Permutations
    {
        public static List<Tuple<T, T>> GetAllPairs<T>(List<T> items)
        {
            var pairs = new List<Tuple<T, T>>();
            var firstItem = items.FirstOrDefault();
            if (firstItem == null || items.Count <= 1)
                return pairs;

            var otherItems = items.Except(new[] { firstItem }).ToList();
            pairs.AddRange(otherItems.Select(item => new Tuple<T, T>(firstItem, item)));
            return pairs.Concat(GetAllPairs(otherItems)).ToList();
        }

        public static List<Tuple<T, T>> GetAllPairs<T>(T[] items)
        {
            return GetAllPairs(new List<T>(items));
        }
    }
}