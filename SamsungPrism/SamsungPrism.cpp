#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <float.h>

using namespace std;

struct Point {
    
    double x, y;     // coordinates
    int cluster;     // no default cluster
    double minDist;  // default infinite distance to nearest cluster

    Point() : x(0.0), y(0.0), cluster(-1), minDist(DBL_MAX) {}
    Point(double x, double y) : x(x), y(y), cluster(-1), minDist(DBL_MAX) {}

    /**
     * Computes the (square) euclidean distance between this point and another
     */
    double distance(Point p) {
        
        return (p.x - x) * (p.x - x) + (p.y - y) * (p.y - y);
    
    }

};

/**
 * Reads in the data.csv file into a vector of points
 * @return vector of points
 *
 */
vector<Point> readcsv() {
    
    vector<Point> points;
    string line;
    ifstream file("C:/Users/adith/Desktop/input_data.csv");

    while (getline(file, line)) {
        
        stringstream lineStream(line);
        string bit;
        
        double x, y;
        getline(lineStream, bit, ',');
        x = stof(bit);
       
        getline(lineStream, bit, '\n');
        y = stof(bit);

        points.push_back(Point(x, y));
    
    }
    
    return points;

}

/**
 * Perform k-means clustering
 * @param points - pointer to vector of points
 * @param epochs - number of k means iterations
 * @param k - the number of initial centroids
 * @return vector of centroids
 */
vector<Point> kMeansClusteringTrain(vector<Point>* points, int epochs, int k) {
    
    int n = points->size();

    // Randomly initialise centroids
    // The index of the centroid within the centroids vector
    // represents the cluster label.
    vector<Point> centroids;
    srand(time(0));
    for (int i = 0; i < k; ++i) {
        
        centroids.push_back(points->at(rand() % n));
    
    }

    for (int i = 0; i < epochs; ++i) {
        
        // For each centroid, compute distance from centroid to each point
        // and update point's cluster if necessary
        for (vector<Point>::iterator c = begin(centroids); c != end(centroids); ++c) {
            
            int clusterId = c - begin(centroids);

            for (vector<Point>::iterator it = points->begin();
                
                it != points->end(); ++it) {
                Point p = *it;
                double dist = c->distance(p);
                
                if (dist < p.minDist) {
                    
                    p.minDist = dist;
                    p.cluster = clusterId;
                
                }
                
                *it = p;
            
            }
        
        }

        // Create vectors to keep track of data needed to compute means
        vector<int> nPoints;
        vector<double> sumX, sumY;
        
        for (int j = 0; j < k; ++j) {
            
            nPoints.push_back(0);
            sumX.push_back(0.0);
            sumY.push_back(0.0);
        
        }

        // Iterate over points to append data to centroids
        for (vector<Point>::iterator it = points->begin(); it != points->end(); ++it) {
            
            int clusterId = it->cluster;
            nPoints[clusterId] += 1;
            sumX[clusterId] += it->x;
            sumY[clusterId] += it->y;

            it->minDist = DBL_MAX;  // reset distance
        
        }
        
        // Compute the new centroids
        for (vector<Point>::iterator c = begin(centroids); c != end(centroids); ++c) {
            
            int clusterId = c - begin(centroids);
            c->x = sumX[clusterId] / nPoints[clusterId];
            c->y = sumY[clusterId] / nPoints[clusterId];
        
        }
    
    }
    
    /*
    // Write to csv
    ofstream myfile;
    myfile.open("C:/Users/adith/Desktop/output.csv");
    myfile << "x,y,c" << endl;

    for (vector<Point>::iterator it = points->begin(); it != points->end(); ++it) {
        
        myfile << it->x << "," << it->y << "," << it->cluster << endl;
    
    }
    
    myfile.close();
    */

    return centroids;
    
}

/**
 * Predicts cluster for a given point
 * @param p - point for which prediction needs to be done
 * @param centroids - pointer to vector of centroids
 * @return cluster id
 */
int kMeansClusteringPredict(Point p, vector<Point>* centroids) {
    
    double min_dist = DBL_MAX;
    int cluster = -1;
    
    //for each centroid find euclidean distance from point
    for (vector<Point>::iterator c = centroids->begin(); c != centroids->end(); ++c) {
        
        int clusterId = c - centroids->begin();
        int clust_dist = c->distance(p);
        cout << "Distance from centroid of cluster " << clusterId << " = " << clust_dist << endl;
        
        //if lower than minimum, then update minimum    
        if (clust_dist < min_dist) {
                
            min_dist = clust_dist;
            cluster = clusterId;
            
        }
    
    }
    
    return cluster;
}

int main() {
    
    vector<Point> points = readcsv();

    // Run k-means with 100 iterations and for 5 clusters
    vector<Point> centroids = kMeansClusteringTrain(&points, 100, 5);

    bool l_end = true;
    Point pred;
    while (l_end) {
        
        char flag;
        
        cout << "Enter point:" << endl;
        cin >> pred.x >> pred.y;
        
        int cluster_id = kMeansClusteringPredict(pred, &centroids);
        cout << "Cluster:" << cluster_id << endl;
        
        cout << "Continue>[y/n]" << endl;
        cin >> flag;
        if ((flag == 'n') || (flag == 'N')) {
            
            l_end = false;
        
        }
    
    }

    cout << "Goodbye!" << endl;

}