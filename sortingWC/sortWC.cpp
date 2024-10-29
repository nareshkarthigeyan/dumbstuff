#include <iostream>
#include <string>
#include <vector>
#include <chrono>
#include <fstream>

using namespace std;

class Algorithm {
    public:
    string name;
    int rank;
    double fastest = 100000;
    double slowest = 0;
    double current = 0;

    Algorithm(string n) {
        name = n;
    }

    void display(){
        cout << "Name: " << name << "\nFastest Time: " << fastest << "\nSlowest Time: " << slowest << "\nCurrent Time: " << current << endl;
    }
};

class BubbleSort : public Algorithm {
    public:
    BubbleSort() : Algorithm("Bubble Sort") {};

    void sort(vector<int> a){
        auto start = chrono::high_resolution_clock::now();

        //Algorithm: 
        int size = a.size();
        for (int i = 0; i < size - 1; i++){
            for(int j = 0; j < size - i - 1; j++){
                if(a[j] > a[j + 1]){
                    int temp = a[j];
                    a[j] = a[j + 1];
                    a[j + 1] = temp;
                }
            }
        }


        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration = end - start;
        double time = duration.count();
        current = time;
        if (time < fastest) fastest = time;
        if (time > slowest) slowest = time;
    }
};

class SelectionSort : public Algorithm {
    public:
    SelectionSort() : Algorithm("Selection Sort") {};

    void sort(vector<int> a){
        auto start = chrono::high_resolution_clock::now();

        //Algorithm: 
        int size = a.size();
        for(int i = 0; i < size - 1; i++){
            int minIndex = i;
            for(int j = i + 1; j < size; j++){
                if(a[j] < a[minIndex]){
                    minIndex = j;
                }
            }
            if(minIndex != i){
                int temp = a[i];
                a[i] = a[minIndex];
                a[minIndex] = temp;
            }
        }

        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration = end - start;
        double time = duration.count();
        current = time;
        if (time < fastest) fastest = time;
        if (time > slowest) slowest = time;
    }
};

class MergeSort : public Algorithm {
    private:
    void merge(vector<int>& arr, int left, int mid, int right){
        int n1 = mid - left + 1;
        int n2 = right - mid;
        vector<int> leftArr(n1);
        vector<int> rightArr(n2);

        for(int i = 0; i < n1; i++) leftArr[i] = arr[left + i];
        for(int j = 0; j < n2; j++) rightArr[j] = arr[mid + 1 + j];

        int i = 0, j = 0, k = left;
        while(i < n1 && j < n2) {
            if(leftArr[i] <= rightArr[j]){
                arr[k] = leftArr[i];
                i++;
            } else {
                arr[j] = rightArr[j];
                j++;
            }
            k++;
        }

        while (i < n1){
            arr[k] = leftArr[i];
            i++;
            k++;
        }

        while (j < n2){
            arr[k] = rightArr[j];
            j++;
            k++;
        }
    }

    void merge_sort(vector<int>& arr, int left, int right) {
        if(left < right){
            int mid = left + (right - left) / 2;

            merge_sort(arr, left, mid);
            merge_sort(arr, mid + 1, right);
            merge(arr, left, mid, right);
        }
    }

    public:
    MergeSort() : Algorithm("Merge Sort") {};

    void sort(vector<int> a){
        auto start = chrono::high_resolution_clock::now();

        //Algorithm: 
        merge_sort(a, 0, a.size() - 1);

        
        auto end = chrono::high_resolution_clock::now();
        chrono::duration<double, milli> duration = end - start;
        double time = duration.count();
        current = time;
        if (time < fastest) fastest = time;
        if (time > slowest) slowest = time;
    }
};

void readFile(const string& filename, vector<int>& arr){
    ifstream file(filename);
    if(!file.is_open()){
        return;
    }

    int num;
    while (file >> num){
        arr.push_back(num);
    }
    file.close();
}

int main(void){
    vector<int> arr;
    cout << arr.size();

    string filename = "input.txt";
    cout << "reading from file:";
    readFile(filename, arr);
    cout << endl << "done" << endl;

    cout << arr.size();

    cout << "bubble sort: " << endl;
    BubbleSort bs;
    bs.sort(arr);
    bs.display();

    cout << "bubble sort: " << endl;
    SelectionSort ss;
    ss.sort(arr);
    ss.display();

    cout << "bubble sort: " << endl;
    MergeSort ms;
    ms.sort(arr);
    ms.display();
}