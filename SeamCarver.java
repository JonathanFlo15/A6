import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.StdOut;

public class SeamCarver {
    //
    private Picture picture;
    //
    private int width;
    //
    private int height;

    // Create a seam carver object based on the given picture
    public SeamCarver(Picture picture) {
        if (picture == null) {
            throw new IllegalArgumentException("Needs valid picture");
        }
        this.picture = new Picture(picture);
        this.width = picture.width();
        this.height = picture.height();
    }

    // Current picture
    public Picture picture() {
        return new Picture(picture);
    }

    // Width of current picture
    public int width() {
        return width;
    }

    // Height of current picture
    public int height() {
        return height;
    }

    // Validate column index
    private void validColumn(int col) {
        if (col < 0 || col >= width)
            throw new IllegalArgumentException("column index out of bounds");
    }

    // Validate row index
    private void validRow(int row) {
        if (row < 0 || row >= height)
            throw new IllegalArgumentException("row index out of bounds");
    }

    // Validate the seam
    private void validSeam(int[] seam, int expectedLength, int validRange) {
        if (seam == null || seam.length != expectedLength)
            throw new IllegalArgumentException("Invalid seam length");
        for (int i = 0; i < seam.length; i++) {
            if (seam[i] < 0 || seam[i] >= validRange)
                throw new IllegalArgumentException("Seam element out of range");
            if (i > 0 && Math.abs(seam[i] - seam[i - 1]) > 1)
                throw new IllegalArgumentException(
                        "Invalid seam (adjacent elements differ by more than 1)");
        }
    }

    private double calculateBorderEnergy(int x, int y) {
        if (x == 0 || x == width - 1) {
            int uy, dy;
            if (y == 0) {
                uy = height - 1;
                dy = y + 1;
            }
            else if (y == height - 1) {
                uy = y - 1;
                dy = 0;
            }
            else {
                uy = y - 1;
                dy = y + 1;
            }

            int rgbUp = picture.getARGB(x, uy);
            int rgbDown = picture.getARGB(x, dy);

            int rUp = (rgbUp >> 16) & 0xFF;
            int gUp = (rgbUp >> 8) & 0xFF;
            int bUp = rgbUp & 0xFF;

            int rDown = (rgbDown >> 16) & 0xFF;
            int gDown = (rgbDown >> 8) & 0xFF;
            int bDown = rgbDown & 0xFF;

            int dyRed = rDown - rUp;
            int dyGreen = gDown - gUp;
            int dyBlue = bDown - bUp;

            return dyRed * dyRed + dyGreen * dyGreen + dyBlue * dyBlue;
        }
        else {
            int lx, rx;
            if (x == 0) {
                lx = width - 1;
                rx = x + 1;
            }
            else {
                lx = x - 1;
                rx = 0;
            }

            int rgbLeft = picture.getARGB(lx, y);
            int rgbRight = picture.getARGB(rx, y);

            int rLeft = (rgbLeft >> 16) & 0xFF;
            int gLeft = (rgbLeft >> 8) & 0xFF;
            int bLeft = rgbLeft & 0xFF;

            int rRight = (rgbRight >> 16) & 0xFF;
            int gRight = (rgbRight >> 8) & 0xFF;
            int bRight = rgbRight & 0xFF;

            int dxRed = rRight - rLeft;
            int dxGreen = gRight - gLeft;
            int dxBlue = bRight - bLeft;

            return dxRed * dxRed + dxGreen * dxGreen + dxBlue * dxBlue;
        }
    }

    public double energy(int x, int y) {
        validColumn(x);
        validRow(y);

        // handle border pixels
        if (x == 0 || x == width - 1 || y == 0 || y == height - 1) {
            return calculateBorderEnergy(x, y);
        }

        // calculate energy for non-border pixels
        int lx = x - 1;
        int rx = x + 1;
        int uy = y - 1;
        int dy = y + 1;

        int rgbLeft = picture.getARGB(lx, y);
        int rgbRight = picture.getARGB(rx, y);
        int rgbUp = picture.getARGB(x, uy);
        int rgbDown = picture.getARGB(x, dy);

        int rLeft = (rgbLeft >> 16) & 0xFF;
        int gLeft = (rgbLeft >> 8) & 0xFF;
        int bLeft = rgbLeft & 0xFF;

        int rRight = (rgbRight >> 16) & 0xFF;
        int gRight = (rgbRight >> 8) & 0xFF;
        int bRight = rgbRight & 0xFF;

        int rUp = (rgbUp >> 16) & 0xFF;
        int gUp = (rgbUp >> 8) & 0xFF;
        int bUp = rgbUp & 0xFF;

        int rDown = (rgbDown >> 16) & 0xFF;
        int gDown = (rgbDown >> 8) & 0xFF;
        int bDown = rgbDown & 0xFF;

        int dxRed = rRight - rLeft;
        int dxGreen = gRight - gLeft;
        int dxBlue = bRight - bLeft;

        int dyRed = rDown - rUp;
        int dyGreen = gDown - gUp;
        int dyBlue = bDown - bUp;

        int dxSquared = dxRed * dxRed + dxGreen * dxGreen + dxBlue * dxBlue;
        int dySquared = dyRed * dyRed + dyGreen * dyGreen + dyBlue * dyBlue;

        return Math.sqrt(dxSquared + dySquared);
    }


    // Energy of pixel at column x and row y
    // public double energy(int x, int y) {
    //     validColumn(x);
    //     validRow(y);
    //
    //
    //     if (x == 0 || x == width - 1 || y == 0 || y == height - 1) {
    //         double borderenergy = calculateBorderEnergy(x, y);
    //         return borderenergy;
    //
    //     }
    //
    //     // Compute energy using the dual-gradient energy function
    //     int leftRGB = picture.getARGB(x - 1, y);
    //     int rightRGB = picture.getARGB(x + 1, y);
    //     int upRGB = picture.getARGB(x, y - 1);
    //     int downRGB = picture.getARGB(x, y + 1);
    //
    //     double deltaX2 = gradient(leftRGB, rightRGB);
    //     double deltaY2 = gradient(upRGB, downRGB);
    //
    //     return Math.sqrt(deltaX2 + deltaY2);
    // }

    // Helper method to compute squared gradient between two RGB values
    private double gradient(int rgb1, int rgb2) {
        int r1 = (rgb1 >> 16) & 0xFF;
        int g1 = (rgb1 >> 8) & 0xFF;
        int b1 = rgb1 & 0xFF;
        int r2 = (rgb2 >> 16) & 0xFF;
        int g2 = (rgb2 >> 8) & 0xFF;
        int b2 = rgb2 & 0xFF;

        return (r2 - r1) * (r2 - r1) +
                (g2 - g1) * (g2 - g1) +
                (b2 - b1) * (b2 - b1);
    }

    // Private method to transpose the image
    private void transpose() {
        Picture transposedPicture = new Picture(height, width);
        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                transposedPicture.setARGB(col, row, picture.getARGB(row, col));
            }
        }
        picture = transposedPicture;

        // Swap width and height
        int temp = width;
        width = height;
        height = temp;
    }

    // Sequence of indices for horizontal seam
    public int[] findHorizontalSeam() {
        transpose();
        int[] hSeam = findVerticalSeam();
        transpose();
        return hSeam;
    }

    // Sequence of indices for vertical seam
    public int[] findVerticalSeam() {
        // Compute energy of each pixel
        double[][] energyMatrix = new double[height][width];
        for (int row = 0; row < height; row++) {
            for (int col = 0; col < width; col++) {
                energyMatrix[row][col] = energy(col, row);
            }
        }

        // Initialize distTo and edgeTo arrays
        double[][] distTo = new double[height][width];
        int[][] edgeTo = new int[height][width];


        // Dynamic programming to compute shortest paths
        for (int col = 0; col < width; col++) {
            distTo[0][col] = energyMatrix[0][col];
        }

        // fill in distTo and edgeTo using dynamic programming
        for (int y = 1; y < height; y++) {
            for (int x = 0; x < width; x++) {
                double minDist = Double.POSITIVE_INFINITY;
                int minCol = -1;

                // check left, center, and right pixels in previous row
                for (int dx = -1; dx <= 1; dx++) {
                    int nx = x + dx;
                    if (nx >= 0 && nx < width) {
                        double dist = distTo[y - 1][nx] + energyMatrix[y][x];
                        if (dist < minDist) {
                            minDist = dist;
                            minCol = nx;
                        }
                    }
                }

                distTo[y][x] = minDist;
                edgeTo[y][x] = minCol;
            }
        }

        // Find the end of the minimum energy seam
        double minTotalEnergy = Double.POSITIVE_INFINITY;
        int minCol = -1;
        for (int col = 0; col < width; col++) {
            if (distTo[height - 1][col] < minTotalEnergy) {
                minTotalEnergy = distTo[height - 1][col];
                minCol = col;
            }
        }

        // Reconstruct the seam from bottom to top
        int[] seam = new int[height];
        seam[height - 1] = minCol;
        for (int row = height - 2; row >= 0; row--) {
            seam[row] = edgeTo[row + 1][seam[row + 1]];
        }

        return seam;
    }


    // Remove horizontal seam from current picture
    public void removeHorizontalSeam(int[] seam) {
        if (height <= 1) {
            throw new IllegalArgumentException("Height must be greater than 1");
        }
        validSeam(seam, width, height);

        transpose();
        removeVerticalSeam(seam);
        transpose();
    }

    // Remove vertical seam from current picture
    public void removeVerticalSeam(int[] seam) {
        if (width <= 1) {
            throw new IllegalArgumentException("Width must be greater than 1");
        }
        validSeam(seam, height, width);

        Picture newPicture = new Picture(width - 1, height);
        for (int row = 0; row < height; row++) {
            int colToRemove = seam[row];
            for (int col = 0; col < colToRemove; col++) {
                newPicture.setARGB(col, row, picture.getARGB(col, row));
            }
            for (int col = colToRemove + 1; col < width - 1; col++) {
                newPicture.setARGB(col, row, picture.getARGB(col + 1, row));
            }
        }
        picture = newPicture;
        width--;
    }

    // Unit testing (optional)
    public static void main(String[] args) {
        Picture picture = new Picture(args[0]);
        SeamCarver seamCarver = new SeamCarver(picture);

        // test width() and height()
        StdOut.println("width: " + seamCarver.width());
        StdOut.println("height: " + seamCarver.height());

        // test energy()
        for (int x = 0; x < seamCarver.width(); x++) {
            for (int y = 0; y < seamCarver.height(); y++) {
                StdOut.printf("energy at (%d, %d): %.2f\n", x, y,
                              seamCarver.energy(x, y));
            }
        }

        // test findVerticalSeam() and removeVerticalSeam()
        int[] verticalSeam = seamCarver.findVerticalSeam();
        StdOut.println("Vertical seam:");
        for (int y : verticalSeam)
            StdOut.print(y + " ");
        StdOut.println();
        seamCarver.removeVerticalSeam(verticalSeam);
        StdOut.println("wdith after vertical seam removal: " + seamCarver.width());

        // test findHorizontalSeam() and removeHorizontalSeam()
        int[] horizontalSeam = seamCarver.findHorizontalSeam();
        StdOut.println("horizonal seam:");
        for (int x : horizontalSeam)
            StdOut.print(x + " ");
        StdOut.println();
        seamCarver.removeHorizontalSeam(horizontalSeam);
        StdOut.println("height after horizontal seam removal: " + seamCarver.height());

        // test picture()
        Picture resultPicture = seamCarver.picture();
        StdOut.println("resulting width: " + resultPicture.width());
        StdOut.println("resulting height: " + resultPicture.height());
    }
}
