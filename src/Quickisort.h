#ifndef _QUICK_I_SORT_H
#define _QUICK_I_SORT_H

/**
 * This script uses the quicksort algorithm for sorting a vector but it
 * only returns the indices of the sorted vector
 */

#include <vector>
#include <utility> // for std::swap

// partition function for quicksort
int partition(std::vector<double>& arr, std::vector<int>& indices, int low, int high, bool descending)
{
    double pivot = arr[high];
    int i = low - 1;
    for (int j = low; j <= high - 1; j++)
    {
        if ((descending && arr[j] >= pivot) || (!descending && arr[j] <= pivot))
        {
            i++;
            std::swap(arr[i], arr[j]);
            std::swap(indices[i], indices[j]);
        }
    }
    std::swap(arr[i + 1], arr[high]);
    std::swap(indices[i + 1], indices[high]);

    return i + 1;
}

// quicksort with index tracking for double type
inline void quicksort(std::vector<double>& arr, std::vector<int>& indices, int low, int high, bool descending)
{
    if (low < high)
    {
        int pi = partition(arr, indices, low, high, descending);
        quicksort(arr, indices, low, pi - 1, descending);
        quicksort(arr, indices, pi + 1, high, descending);
    }
}

// wrapper function to sort array and return sorted indices
inline std::vector<int> qsort_index(const std::vector<double>& arr, bool descending = true)
{
    std::vector<int> indices(arr.size());
    for (size_t i = 0; i < arr.size(); ++i)
    {
        indices[i] = i;
    }
    std::vector<int> sortedIndices = indices; // Copy indices to maintain original order
    quicksort(const_cast<std::vector<double> &>(arr), sortedIndices, 0, arr.size() - 1, descending);

    return sortedIndices;
}

#endif // _QUICK_I_SORT_H
