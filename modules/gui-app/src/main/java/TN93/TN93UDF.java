package TN93;
 
import org.apache.hadoop.hive.ql.exec.UDFArgumentException;
import org.apache.hadoop.hive.ql.exec.UDFArgumentLengthException;
import org.apache.hadoop.hive.ql.metadata.HiveException;
import org.apache.hadoop.hive.ql.udf.generic.GenericUDF;
import org.apache.hadoop.hive.serde2.objectinspector.ObjectInspector;
import org.apache.hadoop.hive.serde2.objectinspector.ObjectInspector.Category;
import org.apache.hadoop.hive.serde2.objectinspector.PrimitiveObjectInspector;
import org.apache.hadoop.hive.serde2.objectinspector.PrimitiveObjectInspector.PrimitiveCategory;
import org.apache.hadoop.hive.serde2.objectinspector.primitive.PrimitiveObjectInspectorConverter;
import org.apache.hadoop.hive.serde2.objectinspector.primitive.PrimitiveObjectInspectorConverter.StringConverter;
import org.apache.hadoop.hive.serde2.objectinspector.primitive.PrimitiveObjectInspectorFactory;
import org.apache.hadoop.hive.serde2.objectinspector.primitive.WritableConstantHiveDecimalObjectInspector;
import org.apache.hadoop.hive.serde2.objectinspector.primitive.WritableConstantStringObjectInspector;
import org.apache.hadoop.io.DoubleWritable;

public class TN93UDF extends GenericUDF {

    private TN93 tn93 = null; 
    private transient StringConverter stringConverter;
    
    /**
     * Arguments are Seq1, Seq1, resolveMode (average, resolve,gapmm,skip), abiguityFraction,
     */
    @Override
    public ObjectInspector initialize(ObjectInspector[] arguments) throws UDFArgumentException {
        if (arguments.length != 4) {
            throw new UDFArgumentLengthException(
                    "TN93() requires 4 arguments (Seq1, Seq1, resolveMode (average, resolve, gapmm, skip), abiguityFraction), got "
                            + arguments.length);
        }
        tn93 = new TN93();

        for (int i = 0; i < arguments.length; i++) {

            if (arguments[i].getCategory() != Category.PRIMITIVE) {
                throw new UDFArgumentException("TN93 only takes primitive types, got " +arguments[i].getCategory() + ":" +  arguments[i].getTypeName());
            }
            PrimitiveObjectInspector primitiveObject = (PrimitiveObjectInspector) arguments[i];
            System.err.println("TN93UDF arg object type : " + primitiveObject.getClass().getName());
        }
        // make sure ambig fraction param is double
        PrimitiveObjectInspector ambigArgOI = (PrimitiveObjectInspector) arguments[3];
        PrimitiveCategory inputType = ambigArgOI.getPrimitiveCategory();

        if (inputType != PrimitiveCategory.DECIMAL) {
            throw new UDFArgumentException("TN93 ambig fraction should be a decimal, got " + arguments[3].getTypeName());
        }
        //initialize TN93 with correct parameters.
        WritableConstantHiveDecimalObjectInspector ambigFractionArg = (WritableConstantHiveDecimalObjectInspector) arguments[3];
        tn93.setMaxAmbiguityFraction(ambigFractionArg.getWritableConstantValue().doubleValue());
        System.err.println("arg3 ambigFractionArg value : " + ambigFractionArg.getWritableConstantValue().doubleValue());
        //TN93UDF arg object type : org.apache.hadoop.hive.serde2.objectinspector.primitive.WritableConstantStringObjectInspector
        WritableConstantStringObjectInspector resolveArg = (WritableConstantStringObjectInspector) arguments[2];
        tn93.setAmbiguityHandling(resolveArg.getWritableConstantValue().toString());
        System.err.println("arg2 resolveArg value : " + resolveArg.getWritableConstantValue().toString());
        
        //string converter for sequence string args.
        PrimitiveObjectInspector sequenceArgsOI = (PrimitiveObjectInspector) arguments[0];
        stringConverter = new PrimitiveObjectInspectorConverter.StringConverter(sequenceArgsOI);
        
        // output object is a double.
        ObjectInspector outputOI = PrimitiveObjectInspectorFactory.writableDoubleObjectInspector;
        return outputOI;
    }

    @Override
    public Object evaluate(DeferredObject[] arguments) throws HiveException {

        String val1 = null;
        String val2 = null;
        if (arguments[0] != null && arguments[1] != null) {
          val1 = (String) stringConverter.convert(arguments[0].get());
          val2 = (String) stringConverter.convert(arguments[1].get());
        }
        if (val1 == null || val2 == null) {
          return null;
        }
        
        Seq seq1 = new Seq("1", val1);
        Seq seq2 = new Seq("2", val2);

        double result = tn93.tn93(seq1,seq2);

        return new DoubleWritable(result);
    }

    @Override
    public String getDisplayString(String[] children) {
        return "TN93 UDF implementation wrapper for SeqRuler.";
    }

}
