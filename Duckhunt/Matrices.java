import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Scanner;
import java.util.stream.Collectors;

public class Matrices {
    static DecimalFormat df = new DecimalFormat("#.000000");
    public static void readMatrix(Scanner scanner, int row, int col, double[][] matrix, boolean display) {
        for (int i = 0; i < row; i++) {
            matrix[i] = new double[col];
            for (int j = 0; j < col; j++) {
                matrix[i][j] = scanner.nextDouble();
                System.out.print(matrix[i][j] + " ");
            }
            System.out.println();
        }
    }

    public static void readMatrix(Scanner scanner, int row, int col, double[][] matrix) {
        for (int i = 0; i < row; i++) {
            matrix[i] = new double[col];
            for (int j = 0; j < col; j++) {
                matrix[i][j] = scanner.nextDouble();
            }
        }
    }

    public static void numByVec(Double num, Double[] vector) {
        for (int i = 0; i < vector.length; i++) {
            vector[i] *= num;
        }
    }
    public static double[] vecByVec(double[] v1, double[] v2) {
        int col1 = v1.length, col2 = v2.length;
        if (col1 != col2) {
            System.err.println("Incompatible vector size!");
            return null;
        }
        double[] result = new double[col1];
        for (int i = 0; i < col1; i++) {
            result[i] = v1[i] * v2[i];
        }
        return result;
    }

    public static double[] vecByMat(double[] vector, double[][] matrix) {
        if (vector == null || matrix == null) {
            System.err.println("Have null params!");
            return null;
        }
        int colV = vector.length, rowM = matrix.length, colM = matrix[0].length;
        if (colV != rowM) {
            try {
                throw new Exception("The dimensions of the vector and matrix do not match");
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        double[] result = new double[colM];
        double sum;
        // Reserve 5 digits after 0
//        java.text.DecimalFormat df = new java.text.DecimalFormat("#.00000");
        for (int j = 0; j < colM; j++) {
            sum = 0.0;
            for (int i = 0; i < rowM; i++) {
                sum += vector[i] * matrix[i][j];
            }
            result[j] = sum;
        }
        return result;
    }

    public static Double elementWise(double[] v1, double[] v2, double[] v3) {
        Double result = 0.;
        for (int i = 0; i < v1.length; i++) {
            result += v1[i] * v2[i] * v3[i];
        }
        return result;
    }

    public static Double[][] transposeMat(Double[][] matrix) {
        int row = matrix.length, col = matrix[0].length;
        Double[][] result = new Double[col][];
        for (int i = 0; i < col; i++) {
            result[i] = new Double[row];
            for (int j = 0; j < row; j++) {
                result[i][j] = matrix[j][i];
            }
        }
        return result;
    }

    public static double[] getMatCol(double[][] matrix, int index) {
        int row = matrix.length, col = matrix[0].length;
        if (index < 0 || index >= col) {
            System.err.println("Illegal column index!");
            return null;
        }
        double[] column = new double[row];
        for (int i = 0; i < row; i++) {
            column[i] = matrix[i][index];
        }
        return column;
    }

    public static double[] getMatRow(double[][] matrix, int index) {
        int row = matrix.length, col = matrix[0].length;
        if (index < 0 || index >= row) {
            System.err.println("Illegal row index!");
            return null;
        }
        double[] result = new double[col];
        for (int i = 0; i < col; i++) {
            result[i] = matrix[index][i];
        }
        return result;
    }

    public static void setMatCol(double[][] matrix, double[] vector, int index) {
        int row = matrix.length;
        if (row != vector.length) {
            System.err.println("Incompatible size!");
            return;
        }
        for (int i = 0; i < row; i++) {
            matrix[i][index] = vector[i];
        }
    }
    public static void setVecZero(Double[] vector) {
        for (int i = 0; i < vector.length; i++) {
            vector[i] = 0.;
        }
    }
    public static void setMatZero(Double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                matrix[i][j] = 0.;
            }
        }
    }

    public static void printVector(double[] vector) {
        if (vector == null) {
            System.err.println("Null vector!");
            return;
        }
        for (int i = 0; i < vector.length; i++) {
            if ((Double)vector[i] < 1E-5) {
                System.err.print(0. + " ");
            } else {
                System.err.print(vector[i] +  "  ");
            }
        }
        System.err.println();
    }

    public static void printMatrix(double[][] matrix) {
        int row = matrix.length, col = matrix[0].length;
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                if (matrix[i][j] < 1E-5) {
                    System.err.print(0. + " ");
                } else{
                    System.err.print(Double.valueOf(df.format(matrix[i][j])) + " ");
                }
            }
            System.err.println();
        }
        System.err.println();
    }

    public static void printCube(double[][][] cube) {
        for (int i = 0; i < cube.length; i++) {
            System.out.println("Dimension t = " + i + ":");
            for (int j = 0; j < cube[0].length; j++) {
                for (int k = 0; k < cube[0][0].length; k++) {
                    System.out.print(cube[i][j][k] + " ");
                }
                System.out.println();
            }
            System.out.println();
        }
    }

    public static void normalize(double[] array) {
        double sum = Arrays.stream(array).sum();
        for (int i = 0; i < array.length; i++) {
            array[i] *= 1.0 / sum;
        }
    }
    public static void normalize(double[][] matrix) {
        for (int i = 0; i < matrix.length; i++) {
            double sum = 0.;
            for (int j = 0; j < matrix[0].length; j++) {
                sum += matrix[i][j];
            }
            sum = 1 / sum;
            for (int j = 0; j < matrix[0].length; j++) {
                matrix[i][j] *= sum;
            }
        }
    }

}
