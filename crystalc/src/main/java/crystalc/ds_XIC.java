/*
 * Licensed to the Apache Software Foundation (ASF) under one
 * or more contributor license agreements.  See the NOTICE file
 * distributed with this work for additional information
 * regarding copyright ownership.  The ASF licenses this file
 * to you under the Apache License, Version 2.0 (the
 * "License"); you may not use this file except in compliance
 * with the License.  You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing,
 * software distributed under the License is distributed on an
 * "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
 * KIND, either express or implied.  See the License for the
 * specific language governing permissions and limitations
 * under the License.
 */

package crystalc;

public class ds_XIC  implements Comparable<ds_XIC>
{
    public double Mz = 0;
    public double Rt = 0;
    public double Height = 0;
    public double Area = 0;
    public double StartRt = 999999;
    public double EndRt = 0;

    public int compareTo(ds_XIC xic)
    {
        return Double.compare(Mz, xic.Mz);
    }
}