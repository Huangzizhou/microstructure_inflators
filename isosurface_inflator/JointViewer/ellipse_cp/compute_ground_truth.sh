for b in {1.0,0.99,0.98,0.9,0.8,0.7,0.5,0.25,0.1,0.05,0.001}; do time ./bisection 1 $b 400; cut -f1-2 -f5 bisection_iterates.txt > exact_dist_$b.txt; done
