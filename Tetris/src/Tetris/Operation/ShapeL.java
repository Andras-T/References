package Tetris.Operation;

/**
 * 'L' alakú tetramino
 */
public class ShapeL extends Tetramino{

    /**
     * Létrehoz egy 'L' alakú tetramino-t
     */
    public ShapeL(){
        super();
        int[][] shape = new int[][]{
                {0, 1, 0, 0},
                {0, 1, 0, 0},
                {0, 1, 1, 0},
                {0, 0, 0, 0}};
        changeShape(shape);
    }

    /**
     * Létrehoz egy másolatot
     * @return másolat
     */
    public ShapeL copy(){
        return new ShapeL();
    }

}
