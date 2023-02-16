package Tetris.Operation;

/**
 * 'Z' alakú tetramino
 */
public class ShapeZ extends Tetramino{

    /**
     * Létrehoz egy 'Z' alakú tetramino-t
     */
    public ShapeZ(){
        super();
        int[][] shape = new int[][]{
                {0, 0, 1, 0},
                {0, 1, 1, 0},
                {0, 1, 0, 0},
                {0, 0, 0, 0}};
        changeShape(shape);
    }

    /**
     * Létrehoz egy másolatot
     * @return másolat
     */
    public ShapeZ copy(){
        return new ShapeZ();
    }

}
