package Tetris.Operation;

/**
 * 'S' alakú tetramino
 */
public class ShapeS extends Tetramino{

    /**
     * Létrehoz egy 'S' alakú tetramino-t
     */
    public ShapeS(){
        super();
        int[][] shape = new int[][]{
                {0, 1, 0, 0},
                {0, 1, 1, 0},
                {0, 0, 1, 0},
                {0, 0, 0, 0}};
        changeShape(shape);
    }

    /**
     * Létrehoz egy másolatot
     * @return másolat
     */
    public ShapeS copy(){
        return new ShapeS();
    }

}
