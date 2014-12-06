using Newtonsoft.Json;
using Newtonsoft.Json.Serialization;

namespace mirsynergy
{
    public static class JsonFileWriter
    {
        public static void WriteToFile<T>(T outputModel, string fileName)
        {
            var json = JsonConvert.SerializeObject(outputModel, new JsonSerializerSettings()
            {
                ContractResolver = new CamelCasePropertyNamesContractResolver(),
                StringEscapeHandling = StringEscapeHandling.EscapeHtml
            });
            System.IO.File.WriteAllText(fileName, json);
        }
    }
}