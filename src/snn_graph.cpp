#include "Rcpp.h"
#include <deque>
#include <algorithm>

/* Adapted from bluster package
 * https://bioconductor.org/packages/release/bioc/html/bluster.html
 */

/* The 'rank' version performs the original SNN clustering described by Xu and Su (2015, Bioinformatics).
 * This defines the weight between two nodes as (k - 0.5 * r), where r is the smallest sum of ranks for any node in both NN-sets.
 * The rank is computed separately in each NN set, with each node being 0-rank in its own set.
 */

// [[Rcpp::export(rng=false)]]
Rcpp::List build_snn_rank(Rcpp::IntegerMatrix neighbors) {
  const size_t k=neighbors.ncol();
  const size_t ncells=neighbors.nrow();

  // Building a host table, identifying the reverse relation from nearest neighbours to cells.
  auto mIt=neighbors.begin();
  std::deque<std::deque<std::pair<size_t, int> > > hosts(ncells);
  for (size_t i=1; i<=k; ++i) {
    for (size_t j=0; j<ncells; ++j, ++mIt) {
      hosts[*mIt - 1].push_back(std::make_pair(i, j)); // Getting to 0-based index, keeping 1-based ranks for now.
    }
  }

  std::deque<int> output_pairs;
  std::deque<double> output_weights;
  std::deque<size_t> current_added;
  std::deque<size_t> current_score(ncells);

  for (size_t j=0; j<ncells; ++j) {
    auto rowtmp=neighbors.row(j);
    auto rtIt=rowtmp.begin();

    int cur_neighbor;
    for (size_t i=0; i<=k; ++i) {
      if (i==0) {
        cur_neighbor=j;
      } else {
        // Adding the actual nearest neighbors for cell 'j'.
        cur_neighbor=*rtIt - 1;
        ++rtIt;

        if (static_cast<size_t>(cur_neighbor) < j) { // avoid duplicates from symmetry in the SNN calculations.
          const size_t& currank=i; // +0, as neighbour 'i' is rank 0 with respect to itself.
          size_t& existing_other=current_score[cur_neighbor];
          if (existing_other==0) {
            existing_other=currank;
            current_added.push_back(cur_neighbor);
          } else if (existing_other > currank) {
            existing_other=currank;
          }
        }
      }

      // Adding the cells connected by shared nearest neighbors, again recording the lowest combined rank per neighbor.
      const auto& hosted=hosts[cur_neighbor];
      for (auto hIt=hosted.begin(); hIt!=hosted.end(); ++hIt) {
        const int& othernode=hIt->second;

        if (static_cast<size_t>(othernode) < j) { // avoid duplicates from symmetry in the SNN calculations.
          size_t currank=hIt->first + i;
          size_t& existing_other=current_score[othernode];
          if (existing_other==0) {
            existing_other=currank;
            current_added.push_back(othernode);
          } else if (existing_other > currank) {
            existing_other=currank;
          }
        }
      }
    }

    for (auto othernode : current_added) {
      // Converting to edges.
      output_pairs.push_back(j + 1);
      output_pairs.push_back(othernode + 1);

      // Ensuring that an edge with a positive weight is always reported.
      size_t& otherscore=current_score[othernode];
      double finalscore = static_cast<double>(k) - 0.5 * static_cast<double>(otherscore);
      output_weights.push_back(std::max(finalscore, 1e-6));

      // Resetting all those added to zero.
      otherscore=0;
    }
    current_added.clear();
  }

  Rcpp::IntegerVector pout(output_pairs.begin(), output_pairs.end());
  Rcpp::NumericVector wout(output_weights.begin(), output_weights.end());
  return Rcpp::List::create(pout, wout);
}

/* The 'number' version performs a much simpler SNN clustering.
 * Here, the weight between two nodes is simply the number of shared nodes in both NN-sets.
 * Each node is also included in its own set, yielding a range of [0, k+1] weights.
 */

// [[Rcpp::export(rng=false)]]
Rcpp::List build_snn_number(Rcpp::IntegerMatrix neighbors) {
  const size_t k=neighbors.ncol();
  const size_t ncells=neighbors.nrow();

  // Building a host table, identifying the reverse relation from nearest neighbours to cells.
  auto mIt=neighbors.begin();
  std::deque<std::deque<size_t> > hosts(ncells);
  for (size_t i=1; i<=k; ++i) {
    for (size_t j=0; j<ncells; ++j, ++mIt) {
      hosts[*mIt - 1].push_back(j); // Getting to 0-based index.
    }
  }

  std::deque<int> output_pairs;
  std::deque<double> output_weights;
  std::deque<size_t> current_added;
  std::deque<size_t> current_score(ncells);

  for (size_t j=0; j<ncells; ++j) {
    auto rowtmp=neighbors.row(j);
    auto rtIt=rowtmp.begin();

    int cur_neighbor;
    for (size_t i=0; i<=k; ++i) {
      if (i==0) {
        cur_neighbor=j;
      } else {
        // Adding the actual nearest neighbors for cell 'j'.
        cur_neighbor=*rtIt - 1;
        ++rtIt;

        if (static_cast<size_t>(cur_neighbor) < j) { // avoid duplicates from symmetry in the SNN calculations.
          size_t& existing_other=current_score[cur_neighbor];
          if (existing_other==0) {
            current_added.push_back(cur_neighbor);
          }
          ++existing_other;
        }
      }

      // Adding the cells connected by shared nearest neighbors, recording the number.
      const auto& hosted=hosts[cur_neighbor];
      for (auto hIt=hosted.begin(); hIt!=hosted.end(); ++hIt) {
        const int& othernode=*hIt;

        if (static_cast<size_t>(othernode) < j) { // avoid duplicates from symmetry in the SNN calculations.
          size_t& existing_other=current_score[othernode];
          if (existing_other==0) {
            current_added.push_back(othernode);
          }
          ++existing_other;
        }
      }
    }

    for (auto othernode : current_added) {
      // Converting to edges.
      output_pairs.push_back(j + 1);
      output_pairs.push_back(othernode + 1);

      // Ensuring that an edge with a positive weight is always reported.
      size_t& otherscore=current_score[othernode];
      output_weights.push_back(std::max(static_cast<double>(otherscore), 1e-6));

      // Resetting all those added to zero.
      otherscore=0;
    }
    current_added.clear();
  }

  Rcpp::IntegerVector pout(output_pairs.begin(), output_pairs.end());
  Rcpp::NumericVector wout(output_weights.begin(), output_weights.end());
  return Rcpp::List::create(pout, wout);
}
